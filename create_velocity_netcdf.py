import math
import os
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path

import h5py
import numpy as np
import scipy
import typer
from netCDF4 import Dataset
from pathos.multiprocessing import ProcessPool
from tqdm import tqdm
from typing_extensions import Annotated

temp_dir = tempfile.mkdtemp()
app = typer.Typer(pretty_exceptions_show_locals=False)

GEOGRIDS_MODEL_FILE = None
POINT_FIELD_RESOLUTION = None
CHUNK_SIZE = None
INVALID_VALUE = None
COLUMNS_TO_USE = ["density", "Vs", "Vp"]


def generate_random_id():
    # Generate a random UUID and return it as a string
    return str(uuid.uuid4())


def apply_rotation_points(points, R):
    return np.dot(points, R.T)


def apply_centering_points(points, center):
    return points - center


def revert_rotation_points(points, R):
    return np.dot(points, R)


def revert_centering_points(points, center):
    return points + center


def get_geodetic(x, y, z):
    rad = 6378137.0  # Equatorial radius
    f = 1.0 / 298.257223563  # Flattening factor
    e2 = 2 * f - f**2  # Square of eccentricity

    lon = np.arctan2(y, x)

    # Initial guess for latitude
    p = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1 - e2))  # Initial latitude approximation

    # Iteratively solve for latitude and altitude
    for _ in range(5):  # Converges quickly in a few iterations
        N = rad / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        alt = p / np.cos(lat) - N
        lat = np.arctan2(z + e2 * N * np.sin(lat), p)

    lat_deg = np.degrees(lat)
    lon_deg = np.degrees(lon)
    return lat_deg, lon_deg, alt


def get_value(lat_longs):
    # Create a temporary directory to store the file
    temp_file = generate_random_id()
    file_path = os.path.join(temp_dir, f"{temp_file}")

    # Write data to the file
    with open(f"{file_path}.in", "w") as file:
        for row in lat_longs:
            lat, lon, alt = row
            file.write(f"{lat}  {lon}  {alt}\n")

    subprocess.run(
        [
            "geomodelgrids_query",
            f"--models={GEOGRIDS_MODEL_FILE}",
            f"--points={file_path}.in",
            f"--output={file_path}.out",
            f"--values={','.join(COLUMNS_TO_USE)}",
        ]
    )
    data = read_lat_lon_file(f"{file_path}.out")
    os.remove(f"{file_path}.in")
    os.remove(f"{file_path}.out")
    return data


def read_lat_lon_file(file_path):
    # Define the data structure with column names
    data = np.genfromtxt(
        file_path,
        skip_header=2,  # Skip the header row
        dtype=float,  # Define as float since all columns are numeric
        names=["x0", "x1", "x2"].extend(COLUMNS_TO_USE),
    )
    return data


def points_to_long_lat(points, center, rotation_matrix):
    output = revert_centering_points(points, center)
    output = revert_rotation_points(output, rotation_matrix)
    lat, long, depth = get_geodetic(output[:, 0], output[:, 1], output[:, 2])
    return np.stack((lat, long, depth)).T


def createNetcdf4ParaviewHandle(
    sname: Path, x: np.ndarray, y: np.ndarray, z: np.ndarray, aName: list[str]
):
    # "create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname.with_name(f"{sname.name}_paraview.nc")
    # print("writing " + fname)
    # Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]
    rootgrp = Dataset(fname, "w", format="NETCDF4")
    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)
    rootgrp.createDimension("w", nz)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    vz = rootgrp.createVariable("w", "f4", ("w",))
    vz[:] = z

    output_arrays = []
    for name in aName:
        output_arrays.append(
            rootgrp.createVariable(name, "f4", ("w", "v", "u"))
        )
    return rootgrp, output_arrays


def createNetcdf4SeisSolHandle(
    sname: Path, x: np.ndarray, y: np.ndarray, z: np.ndarray, aName: list[str]
):
    "create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname.with_name(f"{sname.name}_ASAGI.nc")
    # print("writing " + fname)
    # Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]
    rootgrp = Dataset(fname, "w", format="NETCDF4")
    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)
    rootgrp.createDimension("w", nz)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    vz = rootgrp.createVariable("w", "f4", ("w",))
    vz[:] = z

    ldata4 = [(name, "f4") for name in aName]
    ldata8 = [(name, "f4") for name in aName]  # Does this work?
    mattype4 = np.dtype(ldata4)
    mattype8 = np.dtype(ldata8)
    mat_t = rootgrp.createCompoundType(mattype4, "material")
    mat = rootgrp.createVariable("data", mat_t, ("w", "v", "u"))
    return rootgrp, mat, mattype4, mattype8


# TODO this errors out if it gets an array that has a shape of (1,6)
def get_v_values(i, j, k, xg, yg, zg):
    original_shape = xg.shape
    points = np.stack([xg.flatten(), yg.flatten(), zg.flatten()]).T
    latlongs = points_to_long_lat(points, center, rotation_matrix)
    data = get_value(latlongs)
    values = []
    for i in range(data.shape[1] - 3):
        column = data[:, 3 + i]
        values.append(column.reshape(original_shape))
    return i, j, k, values


def get_clean_values(i, j, k, x_chunk, y_chunk, z_chunk):
    xg, yg, zg = np.meshgrid(x_chunk, y_chunk, z_chunk, indexing="ij")
    _, _, _, values = get_v_values(i, j, k, xg, yg, zg)
    return i, j, k, values


@app.command()
def main(
    meta_file: Annotated[
        Path, typer.Argument(help="Path for a meta file of a mesh")
    ],
    output_filename: Annotated[
        Path, typer.Argument(help="Path to output file")
    ],
    geogrids_model_file: Annotated[
        list[Path], typer.Option(help="Path to GeoGrids model file(s)")
    ] = [
        Path("External/USGS_SFCVM_v21-0_regional.h5"),
    ],
    point_field_resolution: Annotated[
        int, typer.Option(help="Resolution of the output netcdf in m")
    ] = 5000,
    chunk_size: Annotated[
        int,
        typer.Option(
            help="Size of the chunks used to split, because we cant hold the entire matrix in memory"
        ),
    ] = 50,
    cores_to_use: Annotated[
        int,
        typer.Option(help="Number of cores to use while building the netcdf"),
    ] = 8,
    invalid_value: Annotated[
        float, typer.Option(help="Missing value")
    ] = -1.0e20,
    replace_invalid: Annotated[
        bool,
        typer.Option(help="Replace invalid values with closest valid value"),
    ] = True,
    clamp_values: Annotated[
        bool, typer.Option(help="Clamp values of density,Vs and Vp to minimum")
    ] = True,
):
    with h5py.File(meta_file, "r") as f:
        global \
            GEOGRIDS_MODEL_FILE, \
            POINT_FIELD_RESOLUTION, \
            KD_TREE_RESOLUTION, \
            CHUNK_SIZE, \
            REPLACE_INVALID_VALUES, \
            INVALID_VALUE, \
            QueryTree, \
            ValidValues, \
            center, \
            rotation_matrix, \
            COLUMNS_TO_USE

        if len(geogrids_model_file) == 0:
            print("Need to specify atleast a single model file")
            exit(1)

        columns_to_use = COLUMNS_TO_USE
        GEOGRIDS_MODEL_FILE = ",".join(
            [str(file) for file in geogrids_model_file]
        )
        POINT_FIELD_RESOLUTION = point_field_resolution
        CHUNK_SIZE = chunk_size
        INVALID_VALUE = invalid_value
        print(f"Using models: {GEOGRIDS_MODEL_FILE}")

        center = f["center"][:]
        rotation_matrix = f["rotation_matrix"][:]
        bounding_box = f.get("bounding_box")[:]

        min_coords = np.min(bounding_box, axis=0)
        max_coords = np.max(bounding_box, axis=0)

        x = np.linspace(
            min_coords[0],
            max_coords[0],
            int((max_coords[0] - min_coords[0]) / point_field_resolution),
        )
        y = np.linspace(
            min_coords[1],
            max_coords[1],
            int((max_coords[0] - min_coords[0]) / point_field_resolution),
        )
        z = np.linspace(
            min_coords[2],
            max_coords[2],
            int((max_coords[0] - min_coords[0]) / point_field_resolution),
        )

        idx = []
        jdx = []
        kdx = []
        x_chunks = []
        y_chunks = []
        z_chunks = []

        rootgrp, vTd = createNetcdf4ParaviewHandle(
            output_filename, x, y, z, ["rho", "mu", "lambda"]
        )
        rootgrp_ss, seissol_mat, mat_type4, mat_type8 = (
            createNetcdf4SeisSolHandle(
                output_filename, x, y, z, ["rho", "mu", "lambda"]
            )
        )
        for i in range(0, len(x), chunk_size):
            for j in range(0, len(y), chunk_size):
                for k in range(0, len(z), chunk_size):
                    x_chunk = x[i : i + chunk_size]
                    y_chunk = y[j : j + chunk_size]
                    z_chunk = z[k : k + chunk_size]
                    idx.append(i)
                    jdx.append(j)
                    kdx.append(k)
                    x_chunks.append(x_chunk)
                    y_chunks.append(y_chunk)
                    z_chunks.append(z_chunk)
                    # print(f"{vTd[0][i:i + chunk_size, j:j + chunk_size, k:k + chunk_size].shape} {x_chunk.shape} {y_chunk.shape} {z_chunk.shape}")

        processPool = ProcessPool(nodes=cores_to_use)
        results = processPool.imap(
            get_clean_values, idx, jdx, kdx, x_chunks, y_chunks, z_chunks
        )

        final_shape = list(vTd[0].shape)
        final_shape.insert(0, 3)
        temp_arrays = np.zeros(final_shape)

        for result in tqdm(
            results, desc="Generating NetCDF", total=len(x_chunks)
        ):
            i, j, k, values = result
            for column_index, column in enumerate(values):
                column = np.einsum("ijk->kji", column)
                temp_arrays[column_index][
                    k : k + chunk_size, j : j + chunk_size, i : i + chunk_size
                ] = column

        if replace_invalid:
            for column_index in range(len(columns_to_use)):
                print(
                    f"Replacing invalid values with valid values for column {column_index}"
                )
                column_data = temp_arrays[column_index][:]
                invalid_mask = column_data == INVALID_VALUE
                indices = scipy.ndimage.distance_transform_edt(
                    invalid_mask, return_indices=True, return_distances=False
                )
                filled_values = column_data[tuple(indices)]
                temp_arrays[column_index][:, :, :] = filled_values

        if clamp_values:
            temp_arrays[0] = np.clip(  # density
                temp_arrays[0], 2550.0, np.finfo(np.float32).max
            )
            temp_arrays[1] = np.clip(  # vs
                temp_arrays[1], 1950.0, np.finfo(np.float32).max
            )
            temp_arrays[2] = np.clip(  # vp
                temp_arrays[2], 3600.0, np.finfo(np.float32).max
            )

        final_array = np.zeros(final_shape)
        density = temp_arrays[0]
        vs = temp_arrays[1]
        vp = temp_arrays[2]

        final_array[0] = density  # rho
        final_array[1] = np.square(vs) * density  # mu
        final_array[2] = (np.square(vp) - 2 * np.square(vs)) * density  # lambda

        if clamp_values:
            final_array[0] = np.clip(
                final_array[0],
                2550.0,
                np.finfo(np.float32).max,
            )
            final_array[1] = np.clip(
                final_array[1],
                2550.0 * math.pow(1950.0, 2),
                np.finfo(np.float32).max,
            )
            final_array[2] = np.clip(
                final_array[2],
                2550.0 * (3600.0**2 - 2 * (1950**2)),
                np.finfo(np.float32).max,
            )

        for i in range(3):
            vTd[i][:] = final_array[i]

        seissol_array = np.zeros(final_shape[1:], dtype=mat_type4)
        seissol_array["rho"] = final_array[0]
        seissol_array["mu"] = final_array[1]
        seissol_array["lambda"] = final_array[2]
        seissol_mat[:] = seissol_array

        rootgrp.close()
        rootgrp_ss.close()
    shutil.rmtree(temp_dir)


if __name__ == "__main__":
    app()
