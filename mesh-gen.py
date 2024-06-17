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
import gmsh
import sys
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_crs_info
from pyproj import Proj, transform


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


def latlong_to_xy(lat_longs, proj):
    """
    Convert latitude and longitude to x, y coordinates using a projection.

    Parameters:
    lat_longs (numpy array): A 2D numpy array of shape (n, 3) with latitude, longitude and a zero z-value.
    proj (pyproj.Proj): A pyproj projection object.

    Returns:
    numpy array: A 2D numpy array of shape (n, 3) with x, y coordinates and a zero z-value.
    """
    lons = lat_longs[:, 0]
    lats = lat_longs[:, 1]
    xs, ys = proj(lons, lats)
    return np.vstack((xs, ys, np.zeros(len(xs)))).T


def read_csv():
    df = pd.read_csv('Calaverasfaultzone.csv', index_col='Name', encoding='cp1252')
    return list(df.itertuples(index=False, name=None))


def generate_quad(point_tags, do_flip=False):
    a = gmsh.model.geo.addLine(point_tags[0], point_tags[1])
    b = gmsh.model.geo.addLine(point_tags[1], point_tags[2])
    c = gmsh.model.geo.addLine(point_tags[2], point_tags[3])
    d = gmsh.model.geo.addLine(point_tags[3], point_tags[0])
    loop_tag = gmsh.model.geo.addCurveLoop([a, b, c, d])
    if do_flip:
        loop_tag = -loop_tag
    surface_tag = gmsh.model.geo.addPlaneSurface([loop_tag])
    return surface_tag, loop_tag, (a, b, c, d)


def main():
    filtered_records = read_csv()
    print(filtered_records)
    to_generate = filtered_records[18:20]
    num_walls = len(to_generate)
    all_lat_longs = []
    for record in to_generate:
        lat_longs = get_points(record)
        all_lat_longs.append(lat_longs)
    all_long_lats = np.vstack(all_lat_longs)
    bounding_box = generate_extended_bounding_box(all_long_lats, 0.05)

    proj = Proj(proj="utm", zone=10, ellps="WGS84")
    all_xy = latlong_to_xy(all_long_lats, proj)

    min_lat, max_lat, min_lon, max_lon = bounding_box

    # Convert bounding box coordinates to XY coordinates
    bounding_box_lat_longs = np.array([[min_lon, min_lat, 0], [max_lon, max_lat, 0]])
    bounding_box_xy = latlong_to_xy(bounding_box_lat_longs, proj)

    center = get_center(all_xy)
    bounding_box_xy = bounding_box_xy - center
    all_xy = all_xy - center

    BOUNDING_BOX_DEPTH = 10 * 1000  # is in meters
    FAULT_DEPTH = 6 * 1000

    # print("Original Bounding Box (lat/lon):", bounding_box)
    # print("Converted Bounding Box (x/y):", bounding_box_xy)
    # print(all_xy)
    # exit()

    model = gmsh.model
    factory = model.geo
    mesh = model.mesh

    bounding_box_point_tags = []

    gmsh.initialize()
    model.add("Fault")

    bounding_box_point_tags.append(factory.addPoint(bounding_box_xy[0][0], bounding_box_xy[0][1], 0))
    bounding_box_point_tags.append(factory.addPoint(bounding_box_xy[1][0], bounding_box_xy[0][1], 0))
    bounding_box_point_tags.append(factory.addPoint(bounding_box_xy[1][0], bounding_box_xy[1][1], 0))
    bounding_box_point_tags.append(factory.addPoint(bounding_box_xy[0][0], bounding_box_xy[1][1], 0))

    top_quad = generate_quad(bounding_box_point_tags)
    ov2 = gmsh.model.geo.extrude([(2, top_quad[0]),], 0, 0, BOUNDING_BOX_DEPTH)

    gmsh.option.setNumber('Mesh.MeshSizeMin', 10)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 1000)
    # Synchronize and generate mesh
    factory.synchronize()
    gmsh.fltk.run()
    gmsh.finalize()


if __name__ == '__main__':
    main()
