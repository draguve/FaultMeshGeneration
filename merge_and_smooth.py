import csv
import math
import random
import string
from enum import Enum

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from fastkml import kml
from geopy import distance
from scipy.interpolate import make_lsq_spline, make_splprep
from shapely import wkt
from typing_extensions import Annotated

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


def convert_from_multi_string(multi_line_string):
    points = multi_line_string[17:-2].split(",")
    lat_longs = []
    for point in points:
        lat_longs.append(
            [float(x) for x in point.split(" ") if x]
        )  # removes empty string
    lat_longs = np.array(lat_longs)
    return lat_longs


def get_points(multi_line_string):
    # MULTILINESTRING((-121.5036249996688 37.03746799973482 0.0, -121.503775000355 37.03769099972591 0.0))
    points = multi_line_string[17:-2].split(",")  # remove the text
    lat_longs = []
    for point in points:
        lat_longs.append(
            [float(x) for x in point.split(" ") if x]
        )  # removes empty string
    lat_longs = np.array(lat_longs)
    return lat_longs


def get_multi_line_string(array):
    points = list(array)
    formatted = []
    for point in points:
        formatted.append(" ".join(str(x) for x in point))
    return f"MULTILINESTRING(({','.join(formatted)}))"


def get_random(N=10):
    return "".join(random.choices(string.ascii_letters, k=N))


def read_csv(filename):
    df = pd.read_csv(filename, encoding="utf-8")
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
        l = (dx**2 + dy**2) ** 0.5
        return ti + l**alpha

    distance_p1_p2 = distance.distance(
        P1[::-1], P2[::-1]
    ).m  # lat longs are for some reason stored in reverse in the strings
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


def catmull_rom_chain(
    points, distance_btw_points=1000, cat_mull_room_alpha=0.5
):
    point_quadruples = (
        (points[idx + d] for d in range(QUADRUPLE_SIZE))
        for idx in range(num_segments(points))
    )
    all_splines = (
        catmull_rom_spline(*pq, distance_btw_points, alpha=cat_mull_room_alpha)
        for pq in point_quadruples
    )
    return flatten(all_splines)


def generate_catmull_rom(points, distance_btw_points, cat_mull_room_alpha):
    start_point = mirror_point(points[1], points[0])
    end_point = mirror_point(points[-2], points[-1])
    extended_points = np.vstack((start_point, points, end_point))
    chain_points = catmull_rom_chain(
        extended_points[:, 0:2],
        distance_btw_points=distance_btw_points,
        cat_mull_room_alpha=cat_mull_room_alpha,
    )
    # assert len(chain_points) == num_segments(extended_points[:, 0:2]) * NUM_POINTS
    chain_points = np.vstack(chain_points)
    z = np.zeros((len(chain_points), 1), dtype=chain_points.dtype)
    chain_points = np.hstack((chain_points, z))
    return np.vstack((start_point, end_point)), chain_points


def get_distances(points: np.ndarray):
    distances = []
    for i in range(1, points.shape[0]):
        p1 = points[i - 1][0:2]
        p2 = points[i][0:2]
        distance_p1_p2 = distance.distance(p1[::-1], p2[::-1]).m
        distances.append(distance_p1_p2)
    return distances


def generate_bspline(
    pts: np.ndarray, distance_btw_points: float, smoothing: float, k: int
):
    distances = get_distances(pts)
    num_points = math.floor(np.sum(distances) / distance_btw_points)
    x, y = pts[:, 0], pts[:, 1]

    distances.insert(0, 0)
    v = np.cumsum(distances)
    u = v / v[-1]

    spline, u = make_splprep([x, y], u=u, k=k, s=smoothing)
    grid = np.linspace(0, 1, num_points)
    new_points = spline(grid)
    output_array = np.zeros((num_points, 3))
    output_array[:, 0:2] = new_points.T
    return None, output_array


def generate_lsq_spline(pts: np.ndarray, distance_btw_points: float, n=3, k=3):
    distances = get_distances(pts)
    num_points = math.floor(np.sum(distances) / distance_btw_points)
    x, y = pts[:, 0], pts[:, 1]

    distances.insert(0, 0)
    v = np.cumsum(distances)
    u = v / v[-1]

    t_int = np.linspace(0, 1, n + 2)[1:-1]
    t_full = np.r_[np.repeat(u[0], k + 1), t_int, np.repeat(u[-1], k + 1)]

    x_spline = make_lsq_spline(u, x, t_full, k)
    y_spline = make_lsq_spline(u, y, t_full, k)

    u_final = np.linspace(0, 1, num_points)
    x_final = x_spline(u_final)
    y_final = y_spline(u_final)

    print(x_final)
    print(y_final)

    output_array = np.zeros((num_points, 3))
    output_array[:, 0] = x_final
    output_array[:, 1] = y_final

    return None, output_array

    # distances.insert(0, 0)
    # v = np.cumsum(distances)
    # u = v / v[-1]


class SmoothingMode(str, Enum):
    disabled = "disabled"
    catmull = "catmull"
    bspline = "bspline"
    lsq = "lsq"


def main(
    input_filename: Annotated[
        str, typer.Argument(help="Path for a csv file, containing the input")
    ],
    csv_output: Annotated[str, typer.Option(help="Output csv")] = "",
    plot: Annotated[
        bool, typer.Option(help="Show final fault before saving")
    ] = False,
    kml_output: Annotated[str, typer.Option(help="Output kml")] = "",
    smoothing: Annotated[
        SmoothingMode,
        typer.Option(
            case_sensitive=False, help="Mode for smoothing the points"
        ),
    ] = SmoothingMode.bspline,
    resolution: Annotated[
        int,
        typer.Option(
            help="Distance(in m) between the points when using smoothing/Resolution of smoothed output"
        ),
    ] = 1000,
    cat_mull_room_alpha: Annotated[
        float,
        typer.Option(
            help="0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.",
            rich_help_panel="Cat-Mul Rom Options",
        ),
    ] = 0.5,
    remove_close_points: Annotated[
        bool, typer.Option(help="Remove points if they're too close")
    ] = True,
    remove_threshold: Annotated[
        float, typer.Option(help="Threshold for removing points in m")
    ] = 5.0,
    print_distance: Annotated[
        bool, typer.Option(help="Print distance's between the points")
    ] = False,
    force_plot_points: Annotated[
        str, typer.Option(help="Plot extra points")
    ] = "",
    lsq_n: Annotated[
        int,
        typer.Option(
            rich_help_panel="LSQ Options",
        ),
    ] = 3,
    lsq_k: Annotated[
        int,
        typer.Option(
            rich_help_panel="LSQ Options",
        ),
    ] = 3,
    bspline_smooth: Annotated[
        float,
        typer.Option(
            help="0 forces it to go through all the points",
            rich_help_panel="B Spline Options",
        ),
    ] = 0,
    bspline_k: Annotated[
        int,
        typer.Option(
            rich_help_panel="B Spline Options",
        ),
    ] = 3,
    save_image: Annotated[
        str,
        typer.Option(
            help="Location to save the image of the plot",
        ),
    ] = "",
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
        id = "".join([i for i in to_join if not i.isdigit()])
        location = int("".join([i for i in to_join if i.isdigit()]))
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
        id = "-".join([str(id) for id in gen_id])
        name = "-".join(gen_name)
        all_points = np.vstack(all_points)
        all_points = remove_duplicate_points(all_points)
        match smoothing:
            case SmoothingMode.disabled:
                extended_points, curve = None, all_points
            case SmoothingMode.catmull:
                extended_points, curve = generate_catmull_rom(
                    all_points, resolution, cat_mull_room_alpha
                )
                curve = curve[0:-2]
            case SmoothingMode.bspline:
                extended_points, curve = generate_bspline(
                    all_points,
                    distance_btw_points=resolution,
                    smoothing=bspline_smooth,
                    k=bspline_k,
                )
            case SmoothingMode.lsq:
                extended_points, curve = generate_lsq_spline(
                    all_points, distance_btw_points=resolution, n=lsq_n, k=lsq_k
                )

        final_num_points += curve.shape[0]
        outputs.append([id, name, "", all_points, extended_points, curve])
    print(f"Final paths have {final_num_points} points")

    if remove_close_points:
        deleted_points = 0
        for output in outputs:
            indexes = [
                0,
            ]
            last_added = 0
            rom = output[-1]
            for i in range(1, rom.shape[0]):
                p1 = rom[last_added][0:2]
                p2 = rom[i][0:2]
                distance_p1_p2 = distance.distance(p1[::-1], p2[::-1]).m
                if distance_p1_p2 > remove_threshold:
                    indexes.append(i)
                    last_added = i
                else:
                    deleted_points += 1
            indexes = np.array(indexes)
            output[-1] = output[-1][indexes]
        print(f"Deleted {deleted_points} points for being too close")

    # Print Distance between points
    if print_distance:
        for output in outputs:
            rom = output[-1]
            for i in range(1, rom.shape[0]):
                p1 = rom[i - 1][0:2]
                p2 = rom[i][0:2]
                distance_p1_p2 = distance.distance(p1[::-1], p2[::-1]).m
                print(distance_p1_p2)

    if plot:
        for output in outputs:
            plt.plot(
                output[3][:, 0],
                output[3][:, 1],
                c="blue",
                linestyle="-",
                linewidth=0.5,
            )
            if smoothing != SmoothingMode.disabled:
                plt.plot(
                    output[5][:, 0], output[5][:, 1], c="red", linewidth=0.5
                )
                if output[4] is not None:
                    plt.plot(
                        output[4][:, 0],
                        output[4][:, 1],
                        linestyle="none",
                        marker="o",
                        c="green",
                    )

        if force_plot_points != "":
            extra_points = np.array(
                convert_from_multi_string(force_plot_points)
            )
            plt.plot(
                extra_points[:, 0],
                extra_points[:, 1],
                c="pink",
                linestyle="-",
                linewidth=0.5,
            )
            for i, data in enumerate(extra_points):
                plt.text(
                    data[0], data[1], f"Point {i}", fontsize=9, color="blue"
                )

        # legend hack
        plt.plot([], [], "blue", label="Raw input")
        if smoothing != SmoothingMode.disabled:
            plt.plot([], [], "red", label="Smoothed")

        plt.legend(loc="best")
        if save_image != "":
            plt.savefig(save_image)
        else:
            plt.show()

    # Create a KML document
    if kml_output != "":
        k = kml.KML()
        folder = kml.Folder()
        k.append(folder)
        for row in outputs:
            geometry = wkt.loads(get_multi_line_string(row[-1]))
            placemark = kml.Placemark(
                name=row[1],
                geometry=geometry,
            )
            folder.append(placemark)
        with open(kml_output, "w") as f:
            f.write(k.to_string(prettyprint=True))

    # write csv
    if csv_output != "":
        with open(csv_output, "w", newline="") as file:
            csvwriter = csv.writer(file)
            csvwriter.writerow(["ID", "Name", "Geom"])
            for row in outputs:
                csvwriter.writerow(
                    [row[0], row[1], get_multi_line_string(row[-1])]
                )


if __name__ == "__main__":
    typer.run(main)
