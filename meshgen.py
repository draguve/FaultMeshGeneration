import sys
from enum import Enum
from pprint import pprint

import click
import gmsh
import h5py
import meshio
import numpy as np
import pandas as pd
import py3dep
import pyvista as pv
import typer
from scipy.spatial import Delaunay
from skimage.util.shape import view_as_windows
from sklearn.decomposition import PCA
from tqdm import tqdm
from typing_extensions import Annotated
from pathlib import Path


def create_folder_if_valid(path_string):
    # Convert the string to a Path object
    folder_path = Path(path_string)

    try:
        # Check if the path is a valid directory path
        if not folder_path.suffix:  # Folders generally do not have a suffix like ".txt"
            # If the path is valid and doesn't exist, create the folder
            folder_path.mkdir(parents=True, exist_ok=True)
            print(f"Folder created or already exists at: {folder_path}")
        else:
            raise ValueError("The path provided is not a valid directory path.")

    except Exception as e:
        print(f"Error: {e}")


def strided4D(arr, arr2, s):
    return view_as_windows(arr, arr2.shape, step=s)


def stride_conv_strided(arr, arr2, s):
    arr4D = strided4D(arr, arr2, s=s)
    return np.tensordot(arr4D, arr2, axes=((2, 3), (0, 1)))


# from vtkmodules.vtkCommonDataModel import vtkImplicitPolyDataDistance
# from vtkmodules.numpy_interface import dataset_adapter as dsa
def project_to_plane(points, plane_normal, plane_point):
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    projection = points - np.dot(points - plane_point, plane_normal)[:, None] * plane_normal
    return projection


def create_wall(long_lats, rotation_matrix, center, up_diff=1, down_diff=1, step_size=1000.0):
    all_cartesian_top = get_cartesian(long_lats[:, 1], long_lats[:, 0], up_diff * 1000)
    all_cartesian_bottom = get_cartesian(long_lats[:, 1], long_lats[:, 0], -down_diff * 1000)
    all_cartesian_top = apply_rotation_points(all_cartesian_top, rotation_matrix)
    all_cartesian_top = apply_centering_points(all_cartesian_top, center)
    all_cartesian_bottom = apply_rotation_points(all_cartesian_bottom, rotation_matrix)
    all_cartesian_bottom = apply_centering_points(all_cartesian_bottom, center)
    number_of_points, dim = all_cartesian_top.shape

    # height subdivisions
    np.max(np.linalg.norm(all_cartesian_top - all_cartesian_bottom, axis=1, keepdims=True)) / step_size
    num_heights = int(
        np.ceil(np.max(np.linalg.norm(all_cartesian_top - all_cartesian_bottom, axis=1, keepdims=True)) / step_size))
    num_steps_in_between = np.ceil(
        np.linalg.norm(all_cartesian_top[0:-1] - all_cartesian_top[1:], axis=1, keepdims=True) / step_size).astype(int)
    num_steps_in_between = np.squeeze(num_steps_in_between)

    assert (dim == 3)
    # all_walls = pv.MultiBlock()
    all_grids = []
    for i in range(number_of_points - 1):
        points = np.vstack(
            (all_cartesian_top[i], all_cartesian_top[i + 1], all_cartesian_bottom[i + 1], all_cartesian_bottom[i]))
        grid_points = generate_grid(points, num_steps_in_between[i], num_heights)[:, 0:-1, :]
        all_grids.append(grid_points)
    all_points = np.concatenate(all_grids, axis=1)
    connectivity = generate_grid_connectivity(all_points)
    return all_points.reshape(-1, 3), connectivity
    # return pv.PolyData(all_points.reshape(-1, 3), connectivity)


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


def generate_grid(points, num_subdivisions_u, num_subdivisions_v):
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


def convert_from_multi_string(multi_line_string):
    points = multi_line_string[17:-2].split(',')
    lat_longs = []
    for point in points:
        lat_longs.append([float(x) for x in point.split(" ") if x])  # removes empty string
    lat_longs = np.array(lat_longs)
    return lat_longs


def get_points(record):
    # last is long,lats
    multi_line_string = record[-1]
    # MULTILINESTRING((-121.5036249996688 37.03746799973482 0.0, -121.503775000355 37.03769099972591 0.0))
    return convert_from_multi_string(multi_line_string)


def get_cartesian(lat_deg, lon_deg, alt):
    # see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html

    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    rad = np.float64(6378137.0)
    # Radius of the Earth (in meters)
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


def generate_grid_connectivity(verts):
    num_points = verts.shape[0] * verts.shape[1]
    indices = np.array(range(num_points))
    indices = indices.reshape(verts.shape[0], verts.shape[1])
    return generate_connectivity_from_grid_indices(indices)


def generate_connectivity_from_grid_indices(indices):
    first = indices[0:-1, 0:-1].reshape(-1)
    second = indices[1:, 0:-1].reshape(-1)
    third = indices[0:-1, 1:].reshape(-1)
    forth = indices[1:, 1:].reshape(-1)
    dimensions = np.full(((indices.shape[0] - 1) * (indices.shape[1] - 1)), 3, dtype=int)
    first_triangles = np.stack((dimensions, first, second, third)).T
    second_triangles = np.stack((dimensions, second, forth, third)).T
    # verts.reshape(-1, 3)
    return np.vstack((first_triangles, second_triangles)).flatten()


def generate_image_from_dem(dem):
    all_heights = []
    all_lats = []
    all_longs = []
    for label, content in dem.to_pandas().items():
        diffs = content.to_numpy()
        diffs = np.nan_to_num(diffs, nan=0.0, posinf=None, neginf=None)
        lats = np.array(content.index)
        longs = np.zeros(lats.shape)
        longs[:] = label
        all_lats.append(lats)
        all_longs.append(longs)
        all_heights.append(diffs)
    return np.stack(all_lats), np.stack(all_longs), np.stack(all_heights)


def image_to_points(lats, longs, diffs, custom_connectivity=False):
    points = []
    x, y = lats.shape
    for i in range(x):
        verts = get_cartesian(lat_deg=lats[i, :], lon_deg=longs[i, :], alt=diffs[i])  # again this is weird
        points.append(verts)
    verts = np.stack(points)

    if custom_connectivity:
        connectivity = generate_grid_connectivity(verts)
        return verts, connectivity
    return verts, None


def get_edges(topograph_grid_points):
    sizes = topograph_grid_points.shape
    size_x = sizes[0]
    size_y = sizes[1]
    top_sides = np.concatenate((topograph_grid_points[0, :], topograph_grid_points[1:, size_y - 1],
                                np.flip(topograph_grid_points[size_x - 1, :], axis=0)[1:],
                                np.flip(topograph_grid_points[:, 0], axis=0)[1:]))
    return top_sides


def generate_extrusion(lats, longs, extrude_surface_to_depth, rotation_matrix, center, topograph_grid_points,
                       top_topo_connectivity, step_size):
    points = []
    x, y = lats.shape
    for i in range(x):
        verts = get_cartesian(lat_deg=lats[i, :], lon_deg=longs[i, :],
                              alt=-extrude_surface_to_depth * 1000)  # again this is weird
        points.append(verts)
    bottom_grid_points = np.stack(points)
    bottom_grid_points = apply_rotation_points(bottom_grid_points, rotation_matrix)
    bottom_grid_points = apply_centering_points(bottom_grid_points, center)
    topograph_grid_points = apply_rotation_points(topograph_grid_points, rotation_matrix)
    topograph_grid_points = apply_centering_points(topograph_grid_points, center)

    top_topo_num_points = topograph_grid_points.shape[0] * topograph_grid_points.shape[1]
    top_topo_indices = np.array(range(top_topo_num_points))
    top_topo_indices = top_topo_indices.reshape(topograph_grid_points.shape[0], topograph_grid_points.shape[1])

    bottom_topo_num_points = bottom_grid_points.shape[0] * bottom_grid_points.shape[1]
    bottom_topo_indices = top_topo_num_points + np.array(range(bottom_topo_num_points))
    bottom_topo_indices = bottom_topo_indices.reshape(bottom_grid_points.shape[0], bottom_grid_points.shape[1])
    bottom_topo_connectivity = generate_connectivity_from_grid_indices(bottom_topo_indices)

    top_sides = get_edges(topograph_grid_points)
    bottom_sides = get_edges(bottom_grid_points)
    num_heights = int(
        np.ceil(np.max(np.linalg.norm(top_sides - bottom_sides, axis=1, keepdims=True)) / step_size))
    u = np.linspace(0, 1, num_heights + 1)[1:-1]  # remove top and bottom
    side_points = (1 - u[:, np.newaxis, np.newaxis]) * top_sides + u[:, np.newaxis, np.newaxis] * bottom_sides
    side_num_points = side_points.shape[0] * side_points.shape[1]
    side_indices = top_topo_num_points + bottom_topo_num_points + np.array(range(side_num_points))
    side_indices = side_indices.reshape(side_points.shape[0], side_points.shape[1])

    top_sides_indices = get_edges(top_topo_indices)
    bottom_sides_indices = get_edges(bottom_topo_indices)
    side_all_indices = np.vstack((top_sides_indices, side_indices, bottom_sides_indices))
    side_all_indices = np.hstack((side_all_indices, side_all_indices[:, 0, np.newaxis]))
    side_connectivity = generate_connectivity_from_grid_indices(side_all_indices)

    return np.vstack((topograph_grid_points.reshape(-1, 3), bottom_grid_points.reshape(-1, 3),
                      side_points.reshape(-1, 3))), np.concatenate(
        (top_topo_connectivity, bottom_topo_connectivity, side_connectivity))


def read_csv(filename):
    # df = pd.read_csv(filename, index_col='Name', encoding='utf-8')
    df = pd.read_csv(open(filename, errors='replace'), index_col="Name")
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


class TopographySolver(str, Enum):
    vtk = "vtk"
    scipy = "scipy"
    custom = "custom"


class ExtrusionSolver(str, Enum):
    custom = "custom"
    pyvista = "pyvista"


def main(
        input_file: Annotated[str, typer.Argument(help="Path for the input file, containing latitude and longitudes")],

        # Options about fault
        fault_output: Annotated[
            str, typer.Option(help="Fault output filepath (needs to be a folder when splitting to multiple files)",
                              rich_help_panel="Fault Options")] = None,
        fault_height: Annotated[int, typer.Option(help="How high in Km should the fault be above topography",
                                                  rich_help_panel="Fault Options")] = 2,
        fault_depth: Annotated[int, typer.Option(help="How deep in Km should the fault be below topography",
                                                 rich_help_panel="Fault Options")] = 4,
        fault_resolution: Annotated[
            float, typer.Option(help="How big should the triangles in the fault be (in m)",
                                rich_help_panel="Fault Options")] = 50,
        split_to_multiple_files: Annotated[
            bool, typer.Option(help="Split faults into multiple files", rich_help_panel="Fault Options")] = False,

        # Options about topography
        topography_output: Annotated[str, typer.Option(help="Topography output filepath",
                                                       rich_help_panel="Topography Options")] = None,
        just_check_res: Annotated[
            bool, typer.Option(help="Just check all the topography resolutions available for a region",
                               rich_help_panel="Topography Options")] = False,
        topography_resolution: Annotated[
            int, typer.Option(help="Set resolution in m to use", rich_help_panel="Topography Options")] = 30,
        topography_step: Annotated[
            int, typer.Option(help="Stride for topography", rich_help_panel="Topography Options")] = 1,
        surrounding_region: Annotated[float, typer.Option(help="How far in Lat longs to make the bounding box",
                                                          rich_help_panel="Topography Options")] = 0.01,
        topo_solver: Annotated[
            TopographySolver, typer.Option(
                help="What solver to use for point cloud (VTK's crashes on bigger point clouds) (custom does not work with chunked download",
                rich_help_panel="Topography Options")] = TopographySolver.custom,
        compare_solver: Annotated[
            bool, typer.Option(help="Compare generated topography to vtk's delaunay impl",
                               rich_help_panel="Topography Options")] = False,
        compare_topo_resample: Annotated[
            bool, typer.Option(help="Show plot comparing the topography before and after convolution",
                               rich_help_panel="Topography Options")] = False,

        # Options about extrusion
        extrusion_solver: Annotated[
            ExtrusionSolver, typer.Option(
                help="What solver to use for to extrude",
                rich_help_panel="Extrusion Options")] = ExtrusionSolver.custom,
        extrude_surface_to_depth: Annotated[
            float, typer.Option(help="Extrude topography mesh to depth (in Km)",
                                rich_help_panel="Extrusion Options")] = 0.0,

        # Bounding Box Options
        bounding_box_output: Annotated[
            str, typer.Option(help="Where to store the bounding box mesh (does not generate if not specified)",
                              rich_help_panel="Bounding Box Options")] = None,
        bb_distance_from_topography: Annotated[
            int, typer.Option(help="How far away should the bb be from the topography",
                              rich_help_panel="Bounding Box Options")] = 1,
        bb_mesh_size: Annotated[
            float, typer.Option(help="size of the bounding box mesh (in m)",
                                rich_help_panel="Bounding Box Options")] = 500,
        bb_depth_below_fault: Annotated[
            float, typer.Option(help="How deep bounding box be from the topography(in Km)",
                                rich_help_panel="Bounding Box Options")] = 2,
        bb_height_above_topography: Annotated[
            float, typer.Option(help="How high bounding box be above the topography(in Km)",
                                rich_help_panel="Bounding Box Options")] = 0.5,
        plot_bb: Annotated[
            bool, typer.Option(help="Show bounding box in gmsh ui before saving",
                               rich_help_panel="Bounding Box Options")] = False,
        force_bb_location: Annotated[
            bool, typer.Option(help="Force the bb location instead of generating from surrounding region",
                               rich_help_panel="Bounding Box Options")] = False,
        min_bb_lat: Annotated[
            float, typer.Option(help="When forcing bb location min bb lat",
                                rich_help_panel="Bounding Box Options")] = 34.154799999999994,
        max_bb_lat: Annotated[
            float, typer.Option(help="When forcing bb location max bb lat",
                                rich_help_panel="Bounding Box Options")] = 41.208890000000004,
        min_bb_long: Annotated[
            float, typer.Option(help="When forcing bb location min bb long",
                                rich_help_panel="Bounding Box Options")] = -125.14526000000001,
        max_bb_long: Annotated[
            float, typer.Option(help="When forcing bb location max bb long",
                                rich_help_panel="Bounding Box Options")] = -119.00254,
        force_gmsh_bb: Annotated[
            str, typer.Option(help="Force bounding box - ensure points are clockwise order",
                              rich_help_panel="Bounding Box Options")] = "",
        # Misc
        plot: Annotated[bool, typer.Option(help="Show fault and topography mesh",
                                           rich_help_panel="Miscellaneous Options")] = False,
        verbose: Annotated[bool, typer.Option(help="Verbose",
                                              rich_help_panel="Miscellaneous Options")] = False,
        fast_path_disabled: Annotated[
            bool, typer.Option(help="Disable fast path (uses meshio only)",
                               rich_help_panel="Miscellaneous Options")] = False,
        meta_data_output: Annotated[
            str, typer.Option(help="Output file name for h5 file storing center and rotational Matrix",
                              rich_help_panel="Miscellaneous Options")] = None,
):
    if topo_solver != TopographySolver.custom and extrude_surface_to_depth != 0.0 and extrusion_solver == ExtrusionSolver.custom:
        print('Cannot use custom extrusion solver without custom topo solver')
        exit()

    fast_path = False
    if not plot and topo_solver == TopographySolver.custom and not fast_path_disabled:
        if extrude_surface_to_depth > 0.0 and extrusion_solver == ExtrusionSolver.custom:
            print("Using Fast Path")
            fast_path = True
        if extrude_surface_to_depth == 0.0:
            print("Using Fast Path")
            fast_path = True

    if split_to_multiple_files and not fast_path:
        print("Cannot split walls into multiple files without fast path")

    to_generate = read_csv(input_file)
    num_walls = len(to_generate)
    all_lat_longs = []

    if verbose:
        pprint(to_generate)

    for record in to_generate:
        lat_longs = get_points(record)
        if len(lat_longs) == 2:
            print(f"adding extra point for {record[0]}")
            midpoint = (lat_longs[0] + lat_longs[1]) / 2
            lat_longs = np.insert(lat_longs, 1, midpoint, axis=0)
        if len(lat_longs) < 2:
            print(f"Check input file, {record[0]} has only {len(lat_longs)} points")
            exit()
        all_lat_longs.append(lat_longs)

    all_long_lats = np.vstack(all_lat_longs)
    if force_bb_location:
        bounding_box = [min_bb_lat, max_bb_lat, min_bb_long, max_bb_long]
    else:
        bounding_box = generate_extended_bounding_box(all_long_lats, surrounding_region)
    lat_long_bb = bounding_box
    print(bounding_box)
    dep3_bounding_box = get_3dep_bbox(bounding_box[0], bounding_box[1], bounding_box[2], bounding_box[3])
    dem_res = py3dep.check_3dep_availability(dep3_bounding_box)
    if just_check_res:
        print(dem_res)
        exit()

    if f"{topography_resolution}m" not in dem_res or not dem_res[f"{topography_resolution}m"]:
        print(dem_res)
        print(f"Resolution {topography_resolution}m not available")
        exit()

    print(f"Downloading topography for region: {dep3_bounding_box}")
    dem = py3dep.get_dem(dep3_bounding_box, topography_resolution)
    lats, longs, diffs = generate_image_from_dem(dem)

    before = None
    if compare_topo_resample:
        before = get_cartesian(lats.flatten(), longs.flatten(), diffs.flatten())

    # reduce topography if required
    if topography_step > 1:
        kernel = np.ones((topography_step, topography_step)) / (topography_step ** 2)
        lats = stride_conv_strided(lats, kernel, topography_step)
        longs = stride_conv_strided(longs, kernel, topography_step)
        diffs = stride_conv_strided(diffs, kernel, topography_step)

    if compare_topo_resample:
        after = get_cartesian(lats.flatten(), longs.flatten(), diffs.flatten())

        plotter = pv.Plotter()
        plotter.add_mesh(pv.PolyData(before), "red", "wireframe", point_size=5)
        plotter.add_mesh(pv.PolyData(after), "blue", "wireframe", point_size=7)
        plotter.show()

    topograph_grid_points, custom_connectivity = image_to_points(lats, longs, diffs,
                                                                 custom_connectivity=topo_solver == TopographySolver.custom)
    topograph_points = topograph_grid_points.reshape(-1, 3)

    print(f"Num points for topography : {topograph_points.shape}")
    rotational_center = get_center(topograph_points)
    rotation_matrix = get_rotation_matrix_from_direction(rotational_center)
    topograph_points = apply_rotation_points(topograph_points, rotation_matrix)
    center = get_center(topograph_points)
    topograph_points = apply_centering_points(topograph_points, center)
    top_topo_points_only = topograph_points
    topo_surface = None
    if topo_solver == TopographySolver.scipy:
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
    elif topo_solver == TopographySolver.custom:
        if not fast_path:
            topo_surface = pv.PolyData(topograph_points, custom_connectivity)
    else:
        topo_points = pv.PolyData(topograph_points)
        topo_surface = topo_points.delaunay_2d(progress_bar=True)

    if extrude_surface_to_depth > 0:
        if extrusion_solver == ExtrusionSolver.pyvista:
            plane = pv.Plane(
                center=(topo_surface.center[0], topo_surface.center[1], -extrude_surface_to_depth * 1000),
                direction=(0, 0, -1),
                i_size=2000 * 1000,  # code to fix this
                j_size=2000 * 1000,
            )
            topo_surface = topo_surface.extrude_trim((0, 0, -1.0), plane).triangulate()
        elif extrusion_solver == ExtrusionSolver.custom:
            if topo_solver == TopographySolver.custom:
                topograph_points, custom_connectivity = generate_extrusion(lats, longs,
                                                                           extrude_surface_to_depth, rotation_matrix,
                                                                           center,
                                                                           topograph_grid_points, custom_connectivity,
                                                                           1000)  # fault res here for now
                if not fast_path:
                    topo_surface = pv.PolyData(topograph_points, custom_connectivity)

    print(f"Generating faults for {len(to_generate)} sections")
    individual_wall_meshes = []
    all_wall_meshes = pv.MultiBlock()
    # for fast path
    accumulated_num_points = 0
    all_wall_points = []
    all_wall_connectivity = []
    for i in tqdm(range(num_walls), desc="Generating Walls"):
        wall_points, wall_connectivity = create_wall(all_lat_longs[i], rotation_matrix, center, up_diff=fault_height,
                                                     down_diff=fault_depth, step_size=fault_resolution)
        all_wall_points.append(wall_points)
        wall = pv.PolyData(np.array(wall_points, copy=True), np.array(wall_connectivity, copy=True))
        if fast_path:
            wall_triangles = wall_connectivity.reshape(-1, 4)[:, 1:]
            wall_triangles = wall_triangles + accumulated_num_points
            all_wall_connectivity.append(wall_triangles)
            accumulated_num_points = accumulated_num_points + wall_points.shape[0]
            if split_to_multiple_files:
                individual_wall_meshes.append(wall)
        else:
            all_wall_meshes.append(wall)

    bounding_box_points = None
    if bounding_box_output is not None:
        print("Generating bounding box")

        height = np.max(diffs) + bb_height_above_topography * 1000

        box_lats = None
        box_longs = None
        if force_gmsh_bb:
            long_lats = convert_from_multi_string(force_gmsh_bb)
            if len(long_lats) != 4:
                print("force_gmsh_bb needs to have only 4 points")
                exit()
            long_lats = np.array(long_lats)
            box_lats = long_lats[:, 1]
            box_longs = long_lats[:, 0]
            box_points = get_cartesian(box_lats, box_longs, height)
            box_points = apply_rotation_points(box_points, rotation_matrix)
            box_points = apply_centering_points(box_points, center)
            gmsh.initialize()
            gmsh.model.add("Bounding Box")
            p1 = gmsh.model.geo.addPoint(box_points[0][0], box_points[0][1], box_points[0][2], bb_mesh_size)
            p2 = gmsh.model.geo.addPoint(box_points[1][0], box_points[1][1], box_points[1][2], bb_mesh_size)
            p3 = gmsh.model.geo.addPoint(box_points[2][0], box_points[2][1], box_points[2][2], bb_mesh_size)
            p4 = gmsh.model.geo.addPoint(box_points[3][0], box_points[3][1], box_points[3][2], bb_mesh_size)
            l1 = gmsh.model.geo.addLine(p1, p2)
            l2 = gmsh.model.geo.addLine(p2, p3)
            l3 = gmsh.model.geo.addLine(p3, p4)
            l4 = gmsh.model.geo.addLine(p4, p1)
            ll = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            rs = gmsh.model.geo.addSurfaceFilling([ll])
            gmsh.model.geo.extrude([(2, rs)], 0, 0, -(bb_depth_below_fault + fault_depth) * 1000)
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)

        else:
            x_index = np.array(
                [bb_distance_from_topography, bb_distance_from_topography, -bb_distance_from_topography - 1,
                 -bb_distance_from_topography - 1])
            y_index = np.array(
                [bb_distance_from_topography, -bb_distance_from_topography - 1, bb_distance_from_topography,
                 -bb_distance_from_topography - 1])
            box_lats = lats[x_index, y_index]
            box_longs = longs[x_index, y_index]
            box_points = get_cartesian(box_lats, box_longs, height)
            box_points = apply_rotation_points(box_points, rotation_matrix)
            box_points = apply_centering_points(box_points, center)

            gmsh.initialize()
            gmsh.model.add("Bounding Box")
            p1 = gmsh.model.geo.addPoint(box_points[0][0], box_points[0][1], box_points[0][2], bb_mesh_size)
            p2 = gmsh.model.geo.addPoint(box_points[1][0], box_points[1][1], box_points[1][2], bb_mesh_size)
            p3 = gmsh.model.geo.addPoint(box_points[2][0], box_points[2][1], box_points[2][2], bb_mesh_size)
            p4 = gmsh.model.geo.addPoint(box_points[3][0], box_points[3][1], box_points[3][2], bb_mesh_size)

            l1 = gmsh.model.geo.addLine(p1, p2)
            l2 = gmsh.model.geo.addLine(p2, p4)
            l3 = gmsh.model.geo.addLine(p4, p3)
            l4 = gmsh.model.geo.addLine(p3, p1)

            # Create Line Loop and Plane Surface
            ll = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            ps = gmsh.model.geo.addPlaneSurface([ll])

            gmsh.model.geo.extrude([(2, ps)], 0, 0, -(bb_depth_below_fault + fault_depth) * 1000)
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)

        if plot_bb:
            gmsh.fltk.run()
        print(f"Saving bounding box mesh {bounding_box_output}")
        gmsh.write(f"{bounding_box_output}")

        # Get the coordinates of the original and extruded points
        all_points = gmsh.model.getEntities(0)  # Get all points (0 stands for point entity)
        point_coords = []
        for point in all_points:
            coord = gmsh.model.getValue(0, point[1], [])
            point_coords.append(coord)

        # Convert the list of coordinates to a numpy array of shape (8, 3)
        bounding_box_points = np.array(point_coords).reshape(8, 3)
        gmsh.finalize()
        # print(bounding_box_points)
    if fast_path:
        if topography_output is not None:
            print(f"Saving topography : {topography_output}")
            topo_cells = [
                ("triangle", custom_connectivity.reshape(-1, 4)[:, 1:]),
            ]
            topo_mesh = meshio.Mesh(
                topograph_points,
                topo_cells
            )
            topo_mesh.write(topography_output)
        if fault_output is not None:
            if split_to_multiple_files:
                # ensure the location is the name of a folder
                create_folder_if_valid(fault_output)

                for i, wall in tqdm(enumerate(individual_wall_meshes), desc="Saving faults"):
                    pv.save_meshio(f"{fault_output}/{to_generate[i][0]}.stl", wall)
            else:
                print(f"Saving faults : {fault_output}")
                fault_cells = [
                    ("triangle", np.vstack(all_wall_connectivity)),
                ]
                fault_mesh = meshio.Mesh(
                    np.vstack(all_wall_points),
                    fault_cells
                )
                fault_mesh.write(fault_output)
    else:
        if plot:
            plotter = pv.Plotter()
            plotter.add_mesh(topo_surface, "red", "wireframe")
            plotter.add_mesh(all_wall_meshes, "blue", "wireframe")
            if compare_solver and topo_solver != TopographySolver.vtk:
                topo_points2 = pv.PolyData(topograph_points)
                topo_surface2 = topo_points2.delaunay_2d(progress_bar=True)
                plotter.add_mesh(topo_surface2, "green", "wireframe", opacity=0.5)
            elif compare_solver:
                print("Cannot compare vtk")
            if bounding_box_points is not None:
                plotter.add_mesh(pv.PolyData(bounding_box_points), "yellow")
            plotter.show()

        if topography_output is not None:
            print(f"Saving topography : {topography_output}")
            pv.save_meshio(topography_output, topo_surface)
        if fault_output is not None:
            print(f"Saving faults : {fault_output}")
            pv.save_meshio(fault_output, all_wall_meshes.combine())

    print(f"Center: {center}")
    print(f"Rotational Matrix: {rotation_matrix}")
    if meta_data_output is not None:
        print(f"Generating meta file at {meta_data_output}.h5")
        with h5py.File(f'{meta_data_output}.h5', 'w') as hf:
            hf.create_dataset('center', data=center)
            hf.create_dataset('rotation_matrix', data=rotation_matrix)
            hf.create_dataset("bounding_box_latlongs", data=np.array(bounding_box))
            meta = hf.create_group("meta")
            meta.attrs["input_command"] = " ".join(sys.argv)

            all_args = click.get_current_context().params
            for key, item in all_args.items():
                meta.attrs[key] = str(item)

            dt = h5py.special_dtype(vlen=str)
            fault_input = np.array(to_generate)
            dset = hf.create_dataset('fault_input', fault_input.shape, dtype=dt)
            dset[:] = fault_input
            hf.create_dataset("bounding_box_lat_long", data=lat_long_bb)
            hf.create_dataset("topo_points", data=top_topo_points_only)
            hf.create_dataset("all_long_lats", data=all_long_lats)
            hf.create_dataset("fault_points", data=np.vstack(all_wall_points))
            hf.create_dataset("fault_connectivity", data=np.vstack(all_wall_connectivity))

            if bounding_box_points is not None:
                hf.create_dataset("bounding_box", data=bounding_box_points)

            # if bounding_box_output is not None:


if __name__ == "__main__":
    typer.run(main)

