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

def create_wall(long_lats, up_diff=1, down_diff=1, R=6371, connectivity_delta=0):
    all_cartesian_top = get_cartesian(long_lats[:, 1], long_lats[:, 0], R + up_diff)
    all_cartesian_bottom = get_cartesian(long_lats[:, 1], long_lats[:, 0], R - down_diff)
    number_of_points, dim = all_cartesian_top.shape
    assert (dim == 3)
    connectivity = []
    for i in range(number_of_points - 1):
        current_top = i + connectivity_delta
        next_top = i + 1 + connectivity_delta
        current_bottom = number_of_points + i + connectivity_delta
        next_bottom = number_of_points + i + 1 + connectivity_delta
        connectivity.extend([4, current_top, next_top, next_bottom, current_bottom])
    return np.vstack((all_cartesian_top, all_cartesian_bottom)), connectivity


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
    for label, content in dem.to_pandas().items():
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


def read_csv():
    df = pd.read_csv('Calaverasfaultzone.csv', index_col='Name', encoding='cp1252')
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


def apply_centering(mesh, center):
    mesh.points = mesh.points - center
    return mesh


filtered_records = read_csv()
print(filtered_records)
radius = 6371
to_generate = filtered_records[18:19]
all_vertices = []
all_connectivity = []
all_lat_longs = []
for record in to_generate:
    lat_longs = get_points(record)
    all_lat_longs.append(lat_longs)
    new_points, dim = lat_longs.shape
    assert (dim == 3)
    vertices, connectivity = create_wall(lat_longs, up_diff=20, down_diff=250, connectivity_delta=0,
                                         R=radius)
    all_vertices.append(vertices)
    all_connectivity.append(connectivity)
all_long_lats = np.vstack(all_lat_longs)
bounding_box = generate_extended_bounding_box(all_long_lats, 0.01)
dep3_bounding_box = get_3dep_bbox(bounding_box[0], bounding_box[1], bounding_box[2], bounding_box[3])
dem_res = py3dep.check_3dep_availability(dep3_bounding_box)
assert (dem_res["30m"])
res = 30
dem = py3dep.get_dem(dep3_bounding_box, res)
topograph_points = image_to_points(dem, step=10)

rotational_center = get_center(np.vstack(all_vertices))
rotation_matrix = get_rotation_matrix_from_direction(rotational_center)

all_fault_walls = pv.MultiBlock()
topo_points = pv.PolyData(topograph_points)
topo_points = apply_rotation(topo_points, rotation_matrix)
center = get_center(topo_points.points)
topo_points = apply_centering(topo_points, center)
# topo_mesh["elevation"] = np.linalg.norm(topograph_points, axis=1)
# topo_mesh = topo_mesh.delaunay_2d()
# meshes.append(topo_mesh)
for i in range(0, len(all_vertices)):
    fault_wall = pv.PolyData(all_vertices[i], all_connectivity[i])  # subtracting the center here center the mesh
    fault_wall = apply_rotation(fault_wall, rotation_matrix)
    fault_wall = apply_centering(fault_wall, center)
    # mesh["elevation"] = np.linalg.norm(all_vertices[i],axis=1)
    fault_wall = fault_wall.triangulate()
    all_fault_walls.append(fault_wall)  # comment this out to view elevation
# meshes = meshes.combine().clean()

plane = pv.Plane(
    center=(topo_points.center[0], topo_points.center[1], -200),  # need to calculate these automatically
    direction=(0, 0, -1),
    i_size=2000,
    j_size=2000,
)
extruded_mesh = topo_points.delaunay_2d().extrude_trim((0, 0, -1.0),
                                                       plane).triangulate()  # clean this up too many in one go
# plotter = pv.Plotter()
# plotter.show_axes()
# plotter.show_grid()
# plotter.add_mesh(plane, style='wireframe', color='black')
# plotter.add_mesh(extruded_mesh, style="wireframe", color="red")
# plotter.add_mesh(topo_points, show_edges=True)
# plotter.add_mesh(all_fault_walls, show_edges=True, color="yellow")
merged = all_fault_walls.combine().extract_surface().clean()
intersection_meshes = extruded_mesh.intersection(all_fault_walls, split_first=True, split_second=True)
intersection = intersection_meshes[0]
earth_surface_with_wall_points = intersection_meshes[1]
wall_points_with_surface_points = intersection_meshes[2]
mesh1 = earth_surface_with_wall_points.combine().extract_surface().clean()
mesh2 = wall_points_with_surface_points.combine().extract_surface().clean()

enclosed_points = mesh2.select_enclosed_points(mesh1, tolerance=0.1, check_surface=True)
inside_points_mask = enclosed_points.point_data['SelectedPoints'].astype(bool)
# inside_edges = mesh2.extract_points(inside_points_mask, adjacent_cells=False)
inside_points2 = pv.PolyData(mesh2.points[inside_points_mask])
inside_points = pv.PolyData(mesh2.points[inside_points_mask]).delaunay_2d()

# plotter = pv.Plotter()
# plotter.add_mesh(mesh1, color='lightblue', opacity=0.5, label='Mesh 1')
# plotter.add_mesh(inside_points, color='red', label='Filtered Mesh 2')

joined = (mesh1 + inside_points)
# plotter.add_mesh(inside_points2, color='red', style="wireframe")

# plotter.add_mesh(joined, color='green', style="wireframe")

tet = tetgen.TetGen(joined)
tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
grid = tet.grid
# plotter.add_mesh(grid,opacity=0.3,show_edges=True)
cells = grid.cells.reshape(-1, 5)[:, 1:]
cell_center = grid.points[cells].mean(1)
mask = cell_center[:, 2] < -100
cell_ind = mask.nonzero()[0]
subgrid = grid.extract_cells(cell_ind)

# advanced plotting
plotter = pv.Plotter()
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_mesh(joined, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tessellated Mesh ', 'black']])
plotter.show()

# plotter.add_mesh(inside_edges, color='green', label='Filtered Mesh 2')
# plotter.show()