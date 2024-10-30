import h5py
import numpy as np
import os
import tempfile
import uuid
import shutil
import subprocess
from tqdm import tqdm
from netCDF4 import Dataset

meta_file = "outputs/BayModel1_final/meta.h5"
detail_model_file = "External/USGS_SFCVM_v21-0_detailed.h5"
regional_model_file = "External/USGS_SFCVM_v21-0_regional.h5"
output_distance_from_topo = "outputs/VelModelTest/data"
point_field_resolution = 2000  #meter
chunk_size = 50
temp_dir = tempfile.mkdtemp()


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
        ["geomodelgrids_query", f"--models={detail_model_file},{regional_model_file}", f"--points={file_path}.in",
         f"--output={file_path}.out", "--values=Vp,Vs,Qp,Qs"])
    data = read_lat_lon_file(f"{file_path}.out")
    return data


def read_lat_lon_file(file_path):
    # Define the data structure with column names
    data = np.genfromtxt(
        file_path,
        skip_header=2,  # Skip the header row
        dtype=float,  # Define as float since all columns are numeric
        names=["x0", "x1", "x2", "Vp", "Vs", "Qp", "Qs"]
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


with h5py.File(meta_file, "r") as f:
    center = f["center"][:]
    rotation_matrix = f["rotation_matrix"][:]
    bounding_box = f.get("bounding_box")[:]
    min_coords = np.min(bounding_box, axis=0)
    max_coords = np.max(bounding_box, axis=0)

    x = np.linspace(min_coords[0], max_coords[0], int((max_coords[0] - min_coords[0]) / point_field_resolution))
    y = np.linspace(min_coords[1], max_coords[1], int((max_coords[0] - min_coords[0]) / point_field_resolution))
    z = np.linspace(min_coords[2], max_coords[2], int((max_coords[0] - min_coords[0]) / point_field_resolution))

    # latlong = points_to_long_lat(bounding_box, center, rotation_matrix)
    # get_value(latlong)

    rootgrp, vTd = createNetcdf4ParaviewHandle(output_distance_from_topo, z, y, x, "velocity")
    # rootgrp_ss, vTd_ss = createNetcdf4SeisSolHandle(output_distance_from_topo, z, y, x, "velocity")

    total_iterations = (len(x) // chunk_size + 1) * (len(y) // chunk_size + 1) * (len(z) // chunk_size + 1)
    progress_bar = tqdm(total=total_iterations, desc="Generating fault_distance")

    # Loop over the grid in chunks
    for i in range(0, len(x), chunk_size):
        for j in range(0, len(y), chunk_size):
            for k in range(0, len(z), chunk_size):
                x_chunk = x[i:i + chunk_size]
                y_chunk = y[j:j + chunk_size]
                z_chunk = z[k:k + chunk_size]

                xg, yg, zg = np.meshgrid(x_chunk, y_chunk, z_chunk, indexing='ij')
                original_shape = xg.shape
                points = np.stack([xg.flatten(), yg.flatten(), zg.flatten()]).T
                latlongs = points_to_long_lat(points, center, rotation_matrix)
                data = get_value(latlongs)
                values = np.zeros((data.shape[0]))
                for idx, row in enumerate(data):
                    values[idx] = row[3]
                values = values.reshape(original_shape)
                values = np.einsum('ijk->kji', values)

                vTd[k:k + chunk_size, j:j + chunk_size, i:i + chunk_size] = values
                # vTd_ss[k:k + chunk_size, j:j + chunk_size, i:i + chunk_size] = values

                # Update tqdm progress bar
                progress_bar.update(1)

    # Close tqdm progress bar
    progress_bar.close()
    rootgrp.close()
    # rootgrp_ss.close()

shutil.rmtree(temp_dir)
