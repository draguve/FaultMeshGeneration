import pandas as pd
import string
import random
import math
import numpy as np
from fastkml import kml
from shapely import wkt


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


def main():
    input_filename = "output.csv"
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
        outputs.append([id, name, "", all_points])

    # Create a KML document
    k = kml.KML()
    # Create a KML folder
    folder = kml.Folder()
    k.append(folder)
    for row in outputs:
        geometry = wkt.loads(get_multi_line_string(row[-1]))
        placemark = kml.Placemark(name=row[1])
        placemark.geometry = geometry
        folder.append(placemark)
    kml_output = 'output.kml'
    with open(kml_output, 'w') as f:
        f.write(k.to_string(prettyprint=True))


if __name__ == '__main__':
    main()
