import pandas as pd
import string
import random
import math
import numpy as np
from fastkml import kml
from shapely import wkt
import matplotlib.pyplot as plt
import csv
import typer
from typing_extensions import Annotated
from geopy import distance
from merge_and_smooth import get_points, read_csv

from pyproj import Proj, transform


# Function to determine UTM zone
def get_utm_zone(lon):
    return int((lon + 180) // 6) + 1


def main(input_filename="outputs/BayModel1/merged_catmul_500m_a0.5.csv",
         tr_lat=39.45517835297083, tr_lon=-120.27386,
         bl_lat=35.72582999999996, bl_lon=-123.68866751003426):
    data = read_csv(input_filename)
    utm_x = []
    utm_y = []
    for line in data:
        line_lat_longs = get_points(line[-1])
        utm_zones = [get_utm_zone(lon_lat[0]) for lon_lat in line_lat_longs]
        for i, lon_lat in enumerate(line_lat_longs):
            proj = Proj(proj='utm', zone=utm_zones[i], ellps='WGS84', preserve_units=False)
            x, y = proj(lon_lat[0], lon_lat[1])
            utm_x.append(x)
            utm_y.append(y)

    proj = Proj(proj='utm', zone=get_utm_zone(tr_lon), ellps='WGS84', preserve_units=False)
    utm_tr_x, utm_tr_y = proj(tr_lon, tr_lat)
    proj = Proj(proj='utm', zone=get_utm_zone(bl_lon), ellps='WGS84', preserve_units=False)
    utm_bl_x, utm_bl_y = proj(bl_lon, bl_lat)

    fig = plt.figure(frameon=False,figsize=(10,10*((utm_tr_x - utm_bl_x) / (utm_tr_y - utm_bl_y))))
    ax = fig.add_axes([0, 0, 1, 1])
    plt.plot(utm_x, utm_y, marker='o', linestyle='-')
    plt.title('Projected UTM Coordinates')
    plt.xlabel('UTM Easting (m)')
    plt.ylabel('UTM Northing (m)')
    plt.grid()

    plt.xlim(utm_bl_x, utm_tr_x)
    plt.ylim(utm_bl_y, utm_tr_y)
    # plt.gca().set_aspect((utm_tr_x - utm_bl_x) / (utm_tr_y - utm_bl_y))

    # Remove grid, legend, and borders
    plt.axis('off')  # Turn off the axis

    # Save the plot as an SVG file
    plt.savefig('utm_coordinates_plot.svg', format='svg',dpi=300)

    # Show the plot (optional)
    plt.show()


if __name__ == '__main__':
    main()
