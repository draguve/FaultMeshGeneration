from pathlib import Path

import h5py
import numpy as np
import typer
from netCDF4 import Dataset
from pathos.multiprocessing import ProcessPool
from sklearn.neighbors import KDTree
from tqdm import tqdm
from typing_extensions import Annotated

app = typer.Typer(pretty_exceptions_show_locals=False)


def get_clean_values(i, j, k, x_chunk, y_chunk, z_chunk, kd_tree):
    xg, yg, zg = np.meshgrid(x_chunk, y_chunk, z_chunk, indexing="ij")
    pts = np.column_stack((xg.ravel(), yg.ravel(), zg.ravel()))
    dists, _ = kd_tree.query(pts, k=1)
    values = dists.reshape(xg.shape)
    return i, j, k, values


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


@app.command()
def main(
    meta_file: Annotated[
        Path, typer.Argument(help="Path for a meta file of a mesh")
    ],
    output_filename: Annotated[
        Path, typer.Argument(help="Path to output file")
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
):
    with h5py.File(meta_file, "r") as f:
        # POINT_FIELD_RESOLUTION = point_field_resolution
        # CHUNK_SIZE = chunk_size

        center = f["center"][:]
        rotation_matrix = f["rotation_matrix"][:]
        bounding_box = f.get("bounding_box")[:]
        topo_points = f.get("topo_points")[:]
        surface_kd_tree = KDTree(topo_points, leaf_size=30)
        print(f"KD Tree with {topo_points.shape} generated")

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
        trees = []

        rootgrp, vTd = createNetcdf4ParaviewHandle(
            output_filename, x, y, z, ["distance"]
        )
        rootgrp_ss, seissol_mat, mat_type4, mat_type8 = (
            createNetcdf4SeisSolHandle(output_filename, x, y, z, ["distance"])
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
                    trees.append(surface_kd_tree)
                    # print(f"{vTd[0][i:i + chunk_size, j:j + chunk_size, k:k + chunk_size].shape} {x_chunk.shape} {y_chunk.shape} {z_chunk.shape}")

        processPool = ProcessPool(nodes=cores_to_use)
        results = processPool.imap(
            get_clean_values,
            idx,
            jdx,
            kdx,
            x_chunks,
            y_chunks,
            z_chunks,
            trees,
        )

        final_shape = list(vTd[0].shape)
        temp_arrays = np.zeros(final_shape)

        for result in tqdm(
            results, desc="Generating NetCDF", total=len(x_chunks)
        ):
            i, j, k, values = result
            column = np.einsum("ijk->kji", values)
            temp_arrays[
                k : k + chunk_size, j : j + chunk_size, i : i + chunk_size
            ] = column

        vTd[0][:] = temp_arrays[:]

        seissol_array = np.zeros(final_shape, dtype=mat_type4)
        seissol_array["distance"] = temp_arrays[:]
        seissol_mat[:] = seissol_array

        rootgrp.close()
        rootgrp_ss.close()


if __name__ == "__main__":
    app()
