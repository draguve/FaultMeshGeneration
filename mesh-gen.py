import folium
import zipfile
import io
from fastkml import kml
import pandas as pd
from io import StringIO
import csv
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import py3dep
import vtk
import tetgen
import matplotlib.colors as mcolors
from vtkmodules.vtkFiltersCore import vtkCutter
from vtkmodules.vtkCommonDataModel import vtkPlane
import pymeshfix
import typer
import os
from pprint import pprint
import requests_cache
from typing_extensions import Annotated
from scipy.spatial import Delaunay
from sklearn.decomposition import PCA


# from vtkmodules.vtkCommonDataModel import vtkImplicitPolyDataDistance
# from vtkmodules.numpy_interface import dataset_adapter as dsa
def project_to_plane(points, plane_normal, plane_point):
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    projection = points - np.dot(points - plane_point, plane_normal)[:, None] * plane_normal
    return projection


def create_wall(long_lats, rotation_matrix, center, up_diff=1, down_diff=1, step_size=1.0):
    all_cartesian_top = get_cartesian(long_lats[:, 1], long_lats[:, 0], up_diff * 1000)
    all_cartesian_bottom = get_cartesian(long_lats[:, 1], long_lats[:, 0], -down_diff * 1000)
    all_cartesian_top = apply_rotation_points(all_cartesian_top, rotation_matrix)
    all_cartesian_top = apply_centering_points(all_cartesian_top, center)
    all_cartesian_bottom = apply_rotation_points(all_cartesian_bottom, rotation_matrix)
    all_cartesian_bottom = apply_centering_points(all_cartesian_bottom, center)
    number_of_points, dim = all_cartesian_top.shape
    assert (dim == 3)
    all_walls = pv.MultiBlock()
    for i in range(number_of_points - 1):
        points = np.vstack(
            (all_cartesian_top[i], all_cartesian_top[i + 1], all_cartesian_bottom[i + 1], all_cartesian_bottom[i]))
        grid_points = generate_grid(points, step_size)
        point_cloud = pv.PolyData(grid_points.reshape(-1, 3))
        all_walls.append(point_cloud)
    all_walls = all_walls.combine()
    all_walls = all_walls.delaunay_2d()
    return all_walls


def calculate_subdivisions(points, target_resolution=1.0):
    # Calculate the Euclidean distances between adjacent points
    d12 = np.linalg.norm(points[1] - points[0])  # Distance P1 to P2
    d23 = np.linalg.norm(points[2] - points[1])  # Distance P2 to P3
    d34 = np.linalg.norm(points[3] - points[2])  # Distance P3 to P4
    d41 = np.linalg.norm(points[0] - points[3])  # Distance P4 to P1

    # Calculate average lengths of opposite sides
    avg_u = (d12 + d34) / 2
    avg_v = (d23 + d41) / 2

    # Determine number of subdivisions based on target resolution
    num_subdivisions_u = int(avg_u / target_resolution)
    num_subdivisions_v = int(avg_v / target_resolution)

    return num_subdivisions_u, num_subdivisions_v


def generate_grid(points, target_resolution=1.0):
    # Calculate the number of subdivisions
    num_subdivisions_u, num_subdivisions_v = calculate_subdivisions(points, target_resolution)

    # Create a grid of (u, v) values
    u = np.linspace(0, 1, num_subdivisions_u + 1)
    v = np.linspace(0, 1, num_subdivisions_v + 1)
    U, V = np.meshgrid(u, v)

    # Extract the points
    P1, P2, P3, P4 = points

    # Perform bilinear interpolation using NumPy's broadcasting
    grid_points = (1 - U)[:, :, np.newaxis] * (1 - V)[:, :, np.newaxis] * P1 + \
                  U[:, :, np.newaxis] * (1 - V)[:, :, np.newaxis] * P2 + \
                  (1 - U)[:, :, np.newaxis] * V[:, :, np.newaxis] * P4 + \
                  U[:, :, np.newaxis] * V[:, :, np.newaxis] * P3

    return grid_points


def get_center(vertices):
    return np.average(vertices, axis=0)


def get_points(record):
    # last is long,lats
    multi_line_string = record[-1][17:-2]
    # MULTILINESTRING((-121.5036249996688 37.03746799973482 0.0, -121.503775000355 37.03769099972591 0.0))
    points = multi_line_string.split(',')
    lat_longs = []
    for point in points:
        lat_longs.append([float(x) for x in point.split(" ") if x])  # removes empty string
    lat_longs = np.array(lat_longs)
    return lat_longs


def get_cartesian(lat, lon, alt):
    # see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html

    rad = np.float64(6378137.0)  # Radius of the Earth (in meters)
    f = np.float64(1.0 / 298.257223563)  # Flattening factor WGS84 Model
    cosLat = np.cos(lat)
    sinLat = np.sin(lat)
    FF = (1.0 - f) ** 2
    C = 1 / np.sqrt(cosLat ** 2 + FF * sinLat ** 2)
    S = C * FF

    x = (rad * C + alt) * cosLat * np.cos(lon)
    y = (rad * C + alt) * cosLat * np.sin(lon)
    z = (rad * S + alt) * sinLat

    return np.vstack((x, y, z)).T


def get_random_points_nearby(lat_longs, num_random_points=20, max_offset=0.004):
    """
    Generate random points near each lat-long in the input array.

    Parameters:
    lat_longs (numpy array): A 2D numpy array of shape (n, 3) with lat-longs and a zero z-value.
    num_random_points (int): The number of random points to generate near each lat-long.
    max_offset (float): Maximum offset for random points in degrees.

    Returns:
    numpy array: A 2D numpy array with random points near each lat-long.
    """
    random_points = []

    for lat_long in lat_longs:
        lat, lon, z = lat_long
        for _ in range(num_random_points):
            random_lat = lat + np.random.uniform(-max_offset, max_offset)
            random_lon = lon + np.random.uniform(-max_offset, max_offset)
            random_points.append([random_lat, random_lon, z])

    return np.array(random_points)


def generate_bounding_box(lat_longs):
    """
    Generate a bounding box for the given lat-long coordinates.

    Parameters:
    lat_longs (numpy array): A 2D numpy array of shape (n, 3) with lat-longs and a zero z-value.

    Returns:
    tuple: A tuple of (min_lat, max_lat, min_lon, max_lon) representing the bounding box.
    """
    # Extract the latitude and longitude columns
    lats = lat_longs[:, 1]
    lons = lat_longs[:, 0]

    # Calculate the min and max for latitude and longitude
    min_lat = np.min(lats)
    max_lat = np.max(lats)
    min_lon = np.min(lons)
    max_lon = np.max(lons)

    return min_lat, max_lat, min_lon, max_lon


def generate_extended_bounding_box(lat_longs, radius):
    """
    Generate an extended bounding box for the given lat-long coordinates.

    Parameters:
    lat_longs (numpy array): A 2D numpy array of shape (n, 3) with lat-longs and a zero z-value.
    radius (float): The radius by which to extend the bounding box in degrees.

    Returns:
    tuple: A tuple of (min_lat, max_lat, min_lon, max_lon) representing the extended bounding box.
    """
    min_lat, max_lat, min_lon, max_lon = generate_bounding_box(lat_longs)

    # Extend the bounding box by the radius
    min_lat -= radius
    max_lat += radius
    min_lon -= radius
    max_lon += radius

    return min_lat, max_lat, min_lon, max_lon


def generate_random_points_in_bbox(min_lat, max_lat, min_lon, max_lon, num_points):
    """
    Generate random points inside the given bounding box.

    Parameters:
    min_lat (float): The minimum latitude of the bounding box.
    max_lat (float): The maximum latitude of the bounding box.
    min_lon (float): The minimum longitude of the bounding box.
    max_lon (float): The maximum longitude of the bounding box.
    num_points (int): The number of random points to generate.

    Returns:
    numpy array: A 2D numpy array with the generated random points.
    """
    lats = np.random.uniform(min_lat, max_lat, num_points)
    lons = np.random.uniform(min_lon, max_lon, num_points)
    z_vals = np.zeros(num_points)  # Keeping z-value as 0
    random_points = np.column_stack((lons, lats, z_vals))

    return random_points


def get_3dep_bbox(min_lat, max_lat, min_lon, max_lon):
    """
    Check the 3DEP data availability for the given bounding box.

    Parameters:
    min_lat (float): The minimum latitude of the bounding box.
    max_lat (float): The maximum latitude of the bounding box.
    min_lon (float): The minimum longitude of the bounding box.
    max_lon (float): The maximum longitude of the bounding box.

    Returns:
    tuple representing the bounding box required by 3dep package.
    """
    return min_lon, min_lat, max_lon, max_lat


def image_to_points(dep, step=50, scale_factor=1):
    points = []
    index = 0
    for label, content in dep.to_pandas().items():
        index = index + 1
        if index % step != 0:
            continue
        long = label
        diffs = content.to_numpy()
        diffs = np.nan_to_num(diffs, nan=0.0, posinf=None, neginf=None)
        lats = np.array(content.index)
        longs = np.zeros(lats.shape)
        longs[:] = label
        verts = get_cartesian(lat=lats, lon=longs, alt=diffs)  # again this is weird
        points.append(verts[::step])
    verts = np.vstack(points)
    return verts


def read_csv(filename):
    df = pd.read_csv(filename, index_col='Name', encoding='cp1252')
    return list(df.itertuples(index=False, name=None))


def get_rotation_matrix_from_direction(direction):
    normal = direction / np.linalg.norm(direction)  # Normalize the normal vector

    # Define the target normal vector (Z-axis)
    target = np.array([0, 0, 1])

    # Compute the rotation axis (cross product of normal and target)
    axis = np.cross(normal, target)
    axis = axis / np.linalg.norm(axis)  # Normalize the rotation axis

    # Compute the rotation angle (dot product of normal and target)
    angle = np.arccos(np.dot(normal, target))

    # Create the rotation matrix using Rodrigues' rotation formula
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    return R


def get_normal_rotation_matrix(mesh):
    # Compute the normal of the mesh
    normal = mesh.point_normals.mean(axis=0)
    return get_rotation_matrix_from_direction(normal)


def apply_rotation(mesh, R):
    # Apply the rotation matrix to the mesh points
    mesh.points = np.dot(mesh.points, R.T)
    return mesh


def apply_rotation_points(points, R):
    return np.dot(points, R.T)


def apply_centering(mesh, center):
    mesh.points = mesh.points - center
    return mesh


def apply_centering_points(points, center):
    return points - center


def chunk_bounding_box(bounding_box, num_chunks):
    min_lon, min_lat, max_lon, max_lat = bounding_box
    lon_step = (max_lon - min_lon) / num_chunks
    lat_step = (max_lat - min_lat) / num_chunks

    sub_bounding_boxes = []
    for i in range(num_chunks):
        for j in range(num_chunks):
            sub_min_lon = min_lon + i * lon_step
            sub_min_lat = min_lat + j * lat_step
            sub_max_lon = sub_min_lon + lon_step
            sub_max_lat = sub_min_lat + lat_step
            sub_bounding_boxes.append((sub_min_lon, sub_min_lat, sub_max_lon, sub_max_lat))

    return sub_bounding_boxes


def main(
        input_file: Annotated[str, typer.Argument(help="Path for the input file, containing latitude and longitudes")],
        fault_output: Annotated[str, typer.Option(help="Fault output filename")] = "faults.stl",
        topography_output: Annotated[str, typer.Option(help="Topography output filename")] = "topography.stl",
        plot: Annotated[bool, typer.Option(help="Show fault and topography mesh")] = False,
        fault_height: Annotated[int, typer.Option(help="How high in Km should the fault be above topography")] = 10,
        fault_depth: Annotated[int, typer.Option(help="How deep in Km should the fault be below topography")] = 100,
        just_check_res: Annotated[
            bool, typer.Option(help="Just check all the topography resolutions available for a region")] = False,
        topography_resolution: Annotated[int, typer.Option(help="Set resolution in m to use")] = 30,
        verbose: Annotated[bool, typer.Option(help="Verbose")] = False,
        surrounding_region: Annotated[float, typer.Option(help="How far in Lat longs to make the bounding box")] = 0.01,
        topography_step: Annotated[int, typer.Option(help="Stride for topography")] = 1,
        save: Annotated[bool, typer.Option(help="Should you save the output meshes or not")] = True,
        fault_resolution: Annotated[float, typer.Option(help="How big should the triangles in the fault be")] = 1000,
        num_chunks_for_topo: Annotated[int, typer.Option(help="How much to split the topography while loading")] = 1,
        use_scipy_delaunay: Annotated[
            bool, typer.Option(help="Use scipy's delaunay impl, with its PCA impl over VTK's (crashes on bigger point clouds)")] = True,
        comparison_delaunay: Annotated[
            bool, typer.Option(help="If using scipy's delaunay impl will show the comparison in the plot")] = False
):
    filtered_records = read_csv(input_file)
    to_generate = filtered_records[:]
    num_walls = len(to_generate)
    all_lat_longs = []

    if verbose:
        pprint(to_generate)

    for record in to_generate:
        lat_longs = get_points(record)
        all_lat_longs.append(lat_longs)

    all_long_lats = np.vstack(all_lat_longs)
    bounding_box = generate_extended_bounding_box(all_long_lats, surrounding_region)
    dep3_bounding_box = get_3dep_bbox(bounding_box[0], bounding_box[1], bounding_box[2], bounding_box[3])
    dem_res = py3dep.check_3dep_availability(dep3_bounding_box)
    if just_check_res:
        print(dem_res)
        exit()

    if f"{topography_resolution}m" not in dem_res or not dem_res[f"{topography_resolution}m"]:
        print(dem_res)
        print(f"Resolution {topography_resolution}m not available")
        exit()

    topograph_points = None
    if num_chunks_for_topo == 1:
        dem = py3dep.get_dem(dep3_bounding_box, topography_resolution)
        topograph_points = image_to_points(dem, step=topography_step)
    else:
        sub_bounding_boxes = chunk_bounding_box(dep3_bounding_box, num_chunks_for_topo)
        all_topograph_points = []
        for sub_box in tqdm(sub_bounding_boxes, "Getting Topography"):
            dem = py3dep.get_dem(sub_box, topography_resolution)
            tp = image_to_points(dem, step=topography_step)
            all_topograph_points.append(tp)
        topograph_points = np.vstack(all_topograph_points)

    print(f"Num points for topography : {topograph_points.shape}")
    rotational_center = get_center(topograph_points)
    rotation_matrix = get_rotation_matrix_from_direction(rotational_center)
    topograph_points = apply_rotation_points(topograph_points, rotation_matrix)
    center = get_center(topograph_points)
    topograph_points = apply_centering_points(topograph_points, center)
    if use_scipy_delaunay:
        # Perform PCA to reduce to 2D
        print(f"Finding optimal plane to project to")
        pca = PCA(n_components=2)
        points_2d = pca.fit_transform(topograph_points)
        print(f"Calculating triangles based on projection")
        dela = Delaunay(points_2d)
        triangles = dela.simplices
        num_triangles, _ = triangles.shape
        connectivity = np.zeros((num_triangles, 1), dtype=int)
        connectivity[:] = 3
        connectivity = np.hstack((connectivity, triangles)).flatten()
        topo_surface = pv.PolyData(topograph_points, connectivity)
    else:
        topo_points = pv.PolyData(topograph_points)
        topo_surface = topo_points.delaunay_2d(progress_bar=True)

    print(f"Generating faults for {len(to_generate)} sections")
    all_wall_meshes = pv.MultiBlock()
    for i in tqdm(range(num_walls), desc="Generating Walls"):
        wall_multiblock = create_wall(all_lat_longs[i], rotation_matrix, center, up_diff=fault_height,
                                      down_diff=fault_depth, step_size=fault_resolution)
        all_wall_meshes.append(wall_multiblock)

    if plot:
        plotter = pv.Plotter()
        plotter.add_mesh(topo_surface, "red", "wireframe")
        plotter.add_mesh(all_wall_meshes, "blue", "wireframe")
        if comparison_delaunay and use_scipy_delaunay:
            topo_points2 = pv.PolyData(topograph_points)
            topo_surface2 = topo_points2.delaunay_2d(progress_bar=True)
            plotter.add_mesh(topo_surface2, "green", "wireframe", opacity=0.5)
        plotter.show()

    if save:
        pv.save_meshio(topography_output, topo_surface)
        pv.save_meshio(fault_output, all_wall_meshes.combine())


if __name__ == "__main__":
    # main("curved_output.csv", plot=True, save=False, num_chunks_for_topo=4,fault_resolution=500)
    typer.run(main)
