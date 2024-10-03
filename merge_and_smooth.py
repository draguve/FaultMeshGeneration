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

QUADRUPLE_SIZE = 4


def remove_duplicate_points(data):
    # Initialize an empty list to store unique rows
    unique_data = []

    # Iterate through each row in the original data
    for row in data:
        # Check if the row is already in unique_data
        if row.tolist() not in unique_data:
            # If not, add it to unique_data
            unique_data.append(row.tolist())

    # Convert unique_data back to a NumPy array
    unique_data = np.array(unique_data)
    return unique_data


def get_points(multi_line_string):
    # MULTILINESTRING((-121.5036249996688 37.03746799973482 0.0, -121.503775000355 37.03769099972591 0.0))
    points = multi_line_string[17:-2].split(',')  # remove the text
    lat_longs = []
    for point in points:
        lat_longs.append([float(x) for x in point.split(" ") if x])  # removes empty string
    lat_longs = np.array(lat_longs)
    return lat_longs


def get_multi_line_string(array):
    points = list(array)
    formatted = []
    for point in points:
        formatted.append(" ".join(str(x) for x in point))
    return f"MULTILINESTRING(({",".join(formatted)}))"


def get_random(N=10):
    return ''.join(random.choices(string.ascii_letters, k=N))


def read_csv(filename):
    df = pd.read_csv(filename, encoding='utf-8')
    return list(df.itertuples(index=False, name=None))


def mirror_point(p1, p2):
    return 2 * p2 - p1


def num_segments(point_chain):
    return len(point_chain) - (QUADRUPLE_SIZE - 1)


def flatten(list_of_lists):
    return [elem for lst in list_of_lists for elem in lst]


def catmull_rom_spline(P0, P1, P2, P3, distance_btw_points=1000, alpha=0.5):
    def tj(ti, pi, pj):
        xi, yi = pi
        xj, yj = pj
        dx, dy = xj - xi, yj - yi
        l = (dx ** 2 + dy ** 2) ** 0.5
        return ti + l ** alpha

    distance_p1_p2 = distance.distance(P1[::-1],
                                       P2[::-1]).m  # lat longs are for some reason stored in reverse in the strings
    num_points = math.ceil(distance_p1_p2 / distance_btw_points)

    t0 = 0.0
    t1 = tj(t0, P0, P1)
    t2 = tj(t1, P1, P2)
    t3 = tj(t2, P2, P3)
    t = np.linspace(t1, t2, num_points).reshape(num_points, 1)

    A1 = (t1 - t) / (t1 - t0) * P0 + (t - t0) / (t1 - t0) * P1
    A2 = (t2 - t) / (t2 - t1) * P1 + (t - t1) / (t2 - t1) * P2
    A3 = (t3 - t) / (t3 - t2) * P2 + (t - t2) / (t3 - t2) * P3
    B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2
    B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3
    points = (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2
    return points


def catmull_rom_chain(points, distance_btw_points=1000, cat_mull_room_alpha=0.5):
    point_quadruples = ((points[idx + d] for d in range(QUADRUPLE_SIZE)) for idx in range(num_segments(points)))
    all_splines = (catmull_rom_spline(*pq, distance_btw_points, alpha=cat_mull_room_alpha) for pq in point_quadruples)
    return flatten(all_splines)


def generate_catmull_rom(points, distance_btw_points, cat_mull_room_alpha):
    start_point = mirror_point(points[1], points[0])
    end_point = mirror_point(points[-2], points[-1])
    extended_points = np.vstack((start_point, points, end_point))
    chain_points = catmull_rom_chain(extended_points[:, 0:2], distance_btw_points=distance_btw_points,
                                     cat_mull_room_alpha=cat_mull_room_alpha)
    # assert len(chain_points) == num_segments(extended_points[:, 0:2]) * NUM_POINTS
    chain_points = np.vstack(chain_points)
    z = np.zeros((len(chain_points), 1), dtype=chain_points.dtype)
    chain_points = np.hstack((chain_points, z))
    return np.vstack((start_point, end_point)), chain_points


def main(
        input_filename: Annotated[str, typer.Argument(help="Path for a csv file, containing the input")],
        csv_output: Annotated[str, typer.Option(help="Output csv")] = "",
        plot: Annotated[bool, typer.Option(help="Show final fault before saving")] = False,
        kml_output: Annotated[str, typer.Option(help="Output kml")] = "",
        disable_smoothing: Annotated[bool, typer.Option(help="Disables Catmull-Rom smoothing")] = False,
        resolution: Annotated[int, typer.Option(
            help="Distance(in m) between the points when using smoothing/Resolution of smoothed output")] = 1000,
        cat_mull_room_alpha: Annotated[float, typer.Option(help="0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.")] = 0.5
):
    data = read_csv(input_filename)
    lines = {}
    for row in data:
        to_join = row[2]
        if isinstance(to_join, float) and math.isnan(to_join):
            id = get_random()
            lines[id] = {}
            lines[id][0] = row
            continue
        id = ''.join([i for i in to_join if not i.isdigit()])
        location = int(''.join([i for i in to_join if i.isdigit()]))
        if id not in lines:
            lines[id] = {}
        lines[id][location] = row

    outputs = []
    final_num_points = 0
    for name, paths in lines.items():
        order = list(paths.keys())
        order.sort()
        gen_name = []
        gen_id = []
        all_points = []
        for section_key in order:
            gen_id.append(paths[section_key][0])
            gen_name.append(paths[section_key][1])
            multiline = paths[section_key][-1]
            points = get_points(multiline)
            all_points.append(points)
        id = '-'.join([str(id) for id in gen_id])
        name = "-".join(gen_name)
        all_points = np.vstack(all_points)
        all_points = remove_duplicate_points(all_points)
        if disable_smoothing:
            extended_points, catmull = None, all_points
        else:
            extended_points, catmull = generate_catmull_rom(all_points, resolution, cat_mull_room_alpha)
        final_num_points += catmull.shape[0]
        outputs.append([id, name, "", all_points, extended_points, catmull])
    print(f"Final paths have {final_num_points} points")

    if plot:
        for output in outputs:
            plt.plot(output[3][:, 0], output[3][:, 1], c="blue", linestyle="-", label="Raw input", linewidth=0.5)
            if not disable_smoothing:
                plt.plot(output[5][:, 0], output[5][:, 1], c="red", linewidth=0.5, label="Smoothed")
                plt.plot(output[4][:, 0], output[4][:, 1], linestyle="none", marker="o", c="green")
        plt.legend(loc="upper right")
        plt.show()

    # Create a KML document
    if kml_output != "":
        k = kml.KML()
        folder = kml.Folder()
        k.append(folder)
        for row in outputs:
            geometry = wkt.loads(get_multi_line_string(row[-1]))
            placemark = kml.Placemark(name=row[1])
            placemark.geometry = geometry
            folder.append(placemark)
        with open(kml_output, 'w') as f:
            f.write(k.to_string(prettyprint=True))

    # write csv
    if csv_output != "":
        with open(csv_output, 'w', newline="") as file:
            csvwriter = csv.writer(file)
            csvwriter.writerow(["ID", "Name", "Geom"])
            for row in outputs:
                csvwriter.writerow([row[0], row[1], get_multi_line_string(row[-1])])


if __name__ == '__main__':
    typer.run(main)
