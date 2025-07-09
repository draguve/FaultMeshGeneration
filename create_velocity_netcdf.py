import os
import shutil
import subprocess
import tempfile
import uuid

import h5py
import numpy as np
import typer
from netCDF4 import Dataset
from pathos.multiprocessing import ProcessPool
from scipy.spatial import KDTree
from tqdm import tqdm
from typing_extensions import Annotated

temp_dir = tempfile.mkdtemp()

QueryTree = None
ValidValues = None
GEOGRIDS_MODEL_FILE = None
POINT_FIELD_RESOLUTION = None
KD_TREE_RESOLUTION = None
CHUNK_SIZE = None
REPLACE_INVALID_VALUES = None
INVALID_VALUE = None
COLUMN_TO_USE = "Vs"


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
    e2 = 2 * f - f ** 2  # Square of eccentricity

    lon = np.arctan2(y, x)

    # Initial guess for latitude
    p = np.sqrt(x ** 2 + y ** 2)
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
        ["geomodelgrids_query", f"--models={GEOGRIDS_MODEL_FILE}", f"--points={file_path}.in",
         f"--output={file_path}.out", f"--values={COLUMN_TO_USE}"])
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
        names=["x0", "x1", "x2", COLUMN_TO_USE]
    )
    return data


def points_to_long_lat(points, center, rotation_matrix):
    output = revert_centering_points(points, center)
    output = revert_rotation_points(output, rotation_matrix)
    lat, long, depth = get_geodetic(output[:, 0], output[:, 1], output[:, 2])
    return np.stack((lat, long, depth)).T


def createNetcdf4ParaviewHandle(sname, x, y, z, aName):
    #"create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname + "_paraview.nc"
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
    vTd = rootgrp.createVariable(aName, "f4", ("u", "v", "w"))
    return rootgrp, vTd


def createNetcdf4SeisSolHandle(sname, x, y, z, aName):
    "create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname + "_ASAGI.nc"
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
    vTd = rootgrp.createVariable("data", "f4", ("u", "v", "w"))
    return rootgrp, vTd


def get_v_values(i, j, k, xg, yg, zg):
    original_shape = xg.shape
    points = np.stack([xg.flatten(), yg.flatten(), zg.flatten()]).T
    latlongs = points_to_long_lat(points, center, rotation_matrix)
    data = get_value(latlongs)
    values = np.zeros((data.shape[0]))
    for idx, row in enumerate(data):
        values[idx] = row[3]
    values = values.reshape(original_shape)
    return i, j, k, values


def get_clean_values(i, j, k, x_chunk, y_chunk, z_chunk):
    xg, yg, zg = np.meshgrid(x_chunk, y_chunk, z_chunk, indexing='ij')
    _, _, _, values = get_v_values(i, j, k, xg, yg, zg)
    if REPLACE_INVALID_VALUES:
        invalid_mask = values == INVALID_VALUE
        # print(f"Shapes {values.shape} {xg.shape} {invalid_mask.shape}")
        invalid_points = np.column_stack((xg[invalid_mask], yg[invalid_mask], zg[invalid_mask]))
        if len(invalid_points) > 0:
            distances, indices = QueryTree.query(invalid_points)
            values[invalid_mask] = ValidValues[indices]
    values = np.einsum('ijk->kji', values)
    return i, j, k, values


def search_tree(bounding_box):
    min_coords = np.min(bounding_box, axis=0)
    max_coords = np.max(bounding_box, axis=0)
    x = np.linspace(min_coords[0], max_coords[0], int((max_coords[0] - min_coords[0]) / KD_TREE_RESOLUTION))
    y = np.linspace(min_coords[1], max_coords[1], int((max_coords[0] - min_coords[0]) / KD_TREE_RESOLUTION))
    z = np.linspace(min_coords[2], max_coords[2], int((max_coords[0] - min_coords[0]) / KD_TREE_RESOLUTION))

    xg, yg, zg = np.meshgrid(x, y, z, indexing='ij')
    _, _, _, values = get_v_values(0, 0, 0, xg, yg, zg)
    stuff_to_keep = values != INVALID_VALUE
    xg, yg, zg = xg[stuff_to_keep], yg[stuff_to_keep], zg[stuff_to_keep]
    values_to_keep = values[stuff_to_keep]
    print(f"Building Tree of shape {x.shape},{y.shape},{z.shape} with {values_to_keep.shape} values")

    # print(values_to_keep.shape)
    # print(xg.shape, yg.shape, zg.shape)
    points = np.vstack((xg, yg, zg)).T
    tree = KDTree(points)
    return tree, values_to_keep


def main(
        meta_file: Annotated[str, typer.Argument(help="Path for a meta file of a mesh")],
        output_filename: Annotated[str, typer.Argument(help="Path to output file")],
        geogrids_model_file: Annotated[str, typer.Argument(
            help="Path to GeoGrids model file(s)")] = "External/USGS_SFCVM_v21-1_detailed.h5,External/USGS_SFCVM_v21-0_regional.h5",
        point_field_resolution: Annotated[int, typer.Option(help="Resolution of the output netcdf in m")] = 5000,
        kd_tree_resolution: Annotated[
            int, typer.Option(help="Resolution of the kdtree for invalid values in m")] = 20000,
        chunk_size: Annotated[
            int, typer.Option(
                help="Size of the chunks used to split, because we cant hold the entire matrix in memory")] = 50,
        cores_to_use: Annotated[
            int, typer.Option(help="Number of cores to use while building the netcdf")] = 8,
        replace_invalid_values: Annotated[
            bool, typer.Option(help="Replace invalid values from closest valid value")] = True,
        invalid_value: Annotated[
            float, typer.Option(help="Missing value")] = -1.e+20,
        column_to_use: Annotated[
            str, typer.Option(help="What column to get from geogrid model")] = "Vs",
):
    with h5py.File(meta_file, 'r') as f:
        global GEOGRIDS_MODEL_FILE, POINT_FIELD_RESOLUTION, KD_TREE_RESOLUTION, CHUNK_SIZE, REPLACE_INVALID_VALUES, INVALID_VALUE, QueryTree, ValidValues, center, rotation_matrix, COLUMN_TO_USE
        COLUMN_TO_USE = column_to_use
        GEOGRIDS_MODEL_FILE = geogrids_model_file
        POINT_FIELD_RESOLUTION = point_field_resolution
        KD_TREE_RESOLUTION = kd_tree_resolution
        CHUNK_SIZE = chunk_size
        REPLACE_INVALID_VALUES = replace_invalid_values
        INVALID_VALUE = invalid_value

        center = f["center"][:]
        rotation_matrix = f["rotation_matrix"][:]
        bounding_box = f.get("bounding_box")[:]

        if replace_invalid_values:
            QueryTree, ValidValues = search_tree(bounding_box)

        min_coords = np.min(bounding_box, axis=0)
        max_coords = np.max(bounding_box, axis=0)

        x = np.linspace(min_coords[0], max_coords[0], int((max_coords[0] - min_coords[0]) / point_field_resolution))
        y = np.linspace(min_coords[1], max_coords[1], int((max_coords[0] - min_coords[0]) / point_field_resolution))
        z = np.linspace(min_coords[2], max_coords[2], int((max_coords[0] - min_coords[0]) / point_field_resolution))

        idx = []
        jdx = []
        kdx = []
        x_chunks = []
        y_chunks = []
        z_chunks = []

        rootgrp, vTd = createNetcdf4ParaviewHandle(output_filename, z, y, x, column_to_use)
        rootgrp_ss,vTd_ss = createNetcdf4SeisSolHandle(output_filename, z, y, x, column_to_use)
        for i in range(0, len(x), chunk_size):
            for j in range(0, len(y), chunk_size):
                for k in range(0, len(z), chunk_size):
                    x_chunk = x[i:i + chunk_size]
                    y_chunk = y[j:j + chunk_size]
                    z_chunk = z[k:k + chunk_size]
                    idx.append(i)
                    jdx.append(j)
                    kdx.append(k)
                    x_chunks.append(x_chunk)
                    y_chunks.append(y_chunk)
                    z_chunks.append(z_chunk)

        processPool = ProcessPool(nodes=cores_to_use)
        results = processPool.imap(get_clean_values, idx, jdx, kdx, x_chunks, y_chunks, z_chunks)

        for result in tqdm(results, desc="Generating NetCDF", total=len(x_chunks)):
            i, j, k, values = result
            vTd[k:k + chunk_size, j:j + chunk_size, i:i + chunk_size] = values
            vTd_ss[k:k + chunk_size, j:j + chunk_size, i:i + chunk_size] = values

        rootgrp.close()
        rootgrp_ss.close()
    shutil.rmtree(temp_dir)


if __name__ == '__main__':
    typer.run(main)
