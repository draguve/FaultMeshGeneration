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


# from vtkmodules.vtkCommonDataModel import vtkImplicitPolyDataDistance
# from vtkmodules.numpy_interface import dataset_adapter as dsa
def project_to_plane(points, plane_normal, plane_point):
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    projection = points - np.dot(points - plane_point, plane_normal)[:, None] * plane_normal
    return projection


def create_wall(long_lats, up_diff=1, down_diff=1, R=6371):
    all_cartesian_top = get_cartesian(long_lats[:, 1], long_lats[:, 0], R + up_diff)
    all_cartesian_bottom = get_cartesian(long_lats[:, 1], long_lats[:, 0], R - down_diff)
    number_of_points, dim = all_cartesian_top.shape
    assert (dim == 3)
    connectivity = []
    for i in range(number_of_points - 1):
        current_top = i
        next_top = i + 1
        current_bottom = number_of_points + i
        next_bottom = number_of_points + i + 1
        connectivity.extend([4, current_top, next_top, next_bottom, current_bottom])
    return np.vstack((all_cartesian_top, all_cartesian_bottom)), connectivity


def create_wall_with_collision(long_lats, center, rotation_matrix, surface_mesh, extruded_mesh, up_diff=1, down_diff=1,
                               R=6371):
    all_cartesian_top = get_cartesian(long_lats[:, 1], long_lats[:, 0], R + up_diff)
    all_cartesian_bottom = get_cartesian(long_lats[:, 1], long_lats[:, 0], R - down_diff)
    number_of_points, dim = all_cartesian_top.shape
    top_fixed = apply_rotation_points(all_cartesian_top, rotation_matrix)
    top_fixed = apply_centering_points(top_fixed, center)
    bottom_fixed = apply_rotation_points(all_cartesian_bottom, rotation_matrix)
    bottom_fixed = apply_centering_points(bottom_fixed, center)
    bounded_top = np.zeros(all_cartesian_top.shape)
    bounded_bottom = np.zeros(all_cartesian_bottom.shape)
    for i in range(number_of_points):
        points, ind = extruded_mesh.ray_trace(top_fixed[i], bottom_fixed[i])
        points2, ind = extruded_mesh.ray_trace(bottom_fixed[i], top_fixed[i])
        assert points.shape[1] == 3 and points.shape[0] == 1
        assert points2.shape[1] == 3 and points2.shape[0] == 1
        bounded_top[i] = points[0]
        bounded_bottom[i] = points2[0]
    assert (dim == 3)
    wall_mesh = pv.MultiBlock()
    intersection_mesh = pv.MultiBlock()
    for i in range(number_of_points - 1):
        current_top = i
        next_top = i + 1
        current_bottom = number_of_points + i
        next_bottom = number_of_points + i + 1
        quad_points = np.vstack((top_fixed[i:i + 2], bottom_fixed[i:i + 2]))
        # Create a plane from the quad points
        # Here we assume the quad points lie on a plane. We calculate the normal and center of the plane
        quad_center = np.mean(quad_points, axis=0)
        quad_normal = np.cross(quad_points[1] - quad_points[0], quad_points[2] - quad_points[0])
        quad_normal = quad_normal / np.linalg.norm(quad_normal)
        intersecting_with_plane = surface_mesh.slice(quad_normal, quad_center)
        points_intersecting_with_plane = intersecting_with_plane.points

        # Create an orthogonal normal vector (you can choose any vector orthogonal to quad_normal)
        # For simplicity, let's choose one based on an arbitrary vector, e.g., [1, 0, 0]
        # Ensure it's not collinear with quad_normal
        arbitrary_vector = np.array([1, 0, 0])
        if np.allclose(np.dot(quad_normal, arbitrary_vector), 0):
            arbitrary_vector = np.array([0, 1, 0])

        # Ensure it's not collinear with quad_normal
        if np.allclose(np.dot(quad_normal, arbitrary_vector), 0):
            # If it's collinear, adjust the arbitrary vector
            arbitrary_vector = np.array([0, 1, 0])
            if np.allclose(np.dot(quad_normal, arbitrary_vector), 0):
                # If it's still collinear (unlikely), choose another orthogonal vector
                arbitrary_vector = np.array([0, 0, 1])

        orthogonal_vector = np.cross(quad_normal, arbitrary_vector)
        orthogonal_vector = orthogonal_vector / np.linalg.norm(orthogonal_vector)

        # Calculate the orthogonal plane center
        orthogonal_plane_center = quad_center
        projected_points = project_to_plane(points_intersecting_with_plane, orthogonal_vector, orthogonal_plane_center)
        projected_top_points = project_to_plane(top_fixed[i:i + 2], orthogonal_vector, orthogonal_plane_center)

        projected_points_polydata = pv.PolyData(projected_points)
        distance_between_top_point = np.linalg.norm(projected_top_points[1] - projected_top_points[0])
        distance_from_A = np.linalg.norm(projected_points - projected_top_points[0],
                                         axis=1) < distance_between_top_point
        distance_from_B = np.linalg.norm(projected_points - projected_top_points[1],
                                         axis=1) < distance_between_top_point

        between_top_points = points_intersecting_with_plane[distance_from_A & distance_from_B]
        wall_vert = np.vstack(
            (bounded_top[i], between_top_points, bounded_top[i + 1], bounded_bottom[i + 1], bounded_bottom[i]))
        n, _ = wall_vert.shape
        wall_connectivity = [n] + list(range(n))
        wall_mesh.append(pv.PolyData(wall_vert, wall_connectivity))
        intersection_mesh.append(
            pv.PolyData(np.vstack((top_fixed[i], top_fixed[i + 1], bottom_fixed[i + 1], bottom_fixed[i])),
                        [4, 0, 1, 2, 3]))
        # plotter = pv.Plotter()
        # plotter.add_mesh(surface_mesh, 'blue', "wireframe")
        # plotter.add_mesh(pv.PolyData(quad_points, [4, 0, 1, 3, 2]), 'r', 'wireframe')
        # plotter.add_mesh(pv.PolyData(intersecting_with_plane.points), "green")
        # plotter.add_mesh(pv.PolyData(np.vstack((bounded_top[i:i + 2], bounded_bottom[i:i + 2]))), "yellow")
        # plotter.add_mesh(projected_points_polydata, "violet")
        # plotter.add_mesh(pv.PolyData(wall_vert, wall_connectivity), "orange")
        # n, dim = between_top_points.shape
        # if n > 0:
        #     plotter.add_mesh(pv.PolyData(between_top_points), "black")
        # plotter.show()
        # data = is_points_in_quad(quad_points, intersecting_with_surface.points)
        # connectivity.extend([4, current_top, next_top, next_bottom, current_bottom])

    return wall_mesh.combine(), intersection_mesh.combine()

    # for i in range(number_of_points - 1):
    #     points = np.vstack((all_cartesian_top[i:i + 2], all_cartesian_bottom[i:i + 2]))
    #     wall_piece = pv.PolyData(points, [4, 0, 1, 3, 2])
    #     wall_piece = apply_rotation(wall_piece, rotation_matrix)
    #     wall_piece = apply_centering(wall_piece, center)
    #     wall_piece.triangulate().extract_surface().clean()
    #
    #     points, ind = extruded_mesh.ray_trace(top_fixed[i], bottom_fixed[i])
    #     points2, ind = extruded_mesh.ray_trace(bottom_fixed[i], top_fixed[i])
    #     intersection1 = pv.PolyData(points2)
    #     intersection2 = pv.PolyData(points)
    #
    #
    #     # fault_wall = fault_wall.triangulate()
    #     # connectivity.extend([4, current_top, next_top, next_bottom, current_bottom])
    # return wall


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


def get_cartesian(lat=None, lon=None, R=6371, diff=0):  # in km:
    # lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    # R = 6371 # radius of the earth
    R = R + diff
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)
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
        lats = np.array(content.index)
        longs = np.zeros(lats.shape)
        longs[:] = label
        verts = get_cartesian(lat=lats, lon=longs, R=6371 * scale_factor, diff=(diffs / 1000))  # again this is weird
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


def main(
        input_file: str,
        fault_output: str = "faults.stl",
        topography_output: str = "topography.stl",
        show_before_saving: bool = False,
        fault_height: int = 10,
        fault_depth: int = 100,
        just_check_res: bool = False,
        topography_resolution: int = 30,
        verbose: bool = False,
        surrounding_region: float = 0.01,
        topography_step: int = 1
):
    filtered_records = read_csv(input_file)
    radius = 6371
    to_generate = filtered_records[:]
    num_walls = len(to_generate)
    all_lat_longs = []

    print(f"Generating faults for {len(to_generate)} sections")
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

    dem = py3dep.get_dem(dep3_bounding_box, topography_resolution)
    topograph_points = image_to_points(dem, step=topography_step)

    rotational_center = get_center(np.vstack(topograph_points))
    rotation_matrix = get_rotation_matrix_from_direction(rotational_center)

    topo_points = pv.PolyData(topograph_points)
    topo_points = apply_rotation(topo_points, rotation_matrix)
    center = get_center(topo_points.points)
    topo_points = apply_centering(topo_points, center)

    topo_surface = topo_points.delaunay_2d()

    all_wall_meshes = pv.MultiBlock()
    for i in range(num_walls):
        wall_verts, wall_connectivity, = create_wall(all_lat_longs[i], up_diff=fault_height, down_diff=fault_depth,
                                                     R=radius)
        wall_verts = apply_rotation_points(wall_verts, rotation_matrix)
        wall_verts = apply_centering_points(wall_verts, center)
        wall_mesh = pv.PolyData(wall_verts, wall_connectivity).triangulate()
        all_wall_meshes.append(wall_mesh)

    if show_before_saving:
        plotter = pv.Plotter()
        plotter.add_mesh(topo_surface, "red", "wireframe")
        plotter.add_mesh(all_wall_meshes, "blue", "wireframe")
        plotter.show()

    pv.save_meshio(topography_output, topo_surface)
    pv.save_meshio(fault_output, all_wall_meshes.combine())


if __name__ == "__main__":
    # main("Calaverasfaultzone.csv", show_before_saving=True)
    typer.run(main)
