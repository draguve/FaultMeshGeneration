import csv

rows = []
with open('raw.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        rows.append(row)

output_rows = []
for i, row in enumerate(rows):
    if i == 0:
        continue  # header
    if row[0] == "":
        continue  # empty

    cleaned_data = [item for item in row[18:] if item]
    multiline_string = "MULTILINESTRING((" + ",".join(
        [f"{lon} {lat} 0.0" for lat, lon in zip(cleaned_data[::2], cleaned_data[1::2])]) + "))"
    id = row[0]
    name = row[1]
    output_rows.append([id, name, "", multiline_string])

with open("ouput.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["id", "Name", "to_join", "multilinestring"])
    writer.writerows(output_rows)
