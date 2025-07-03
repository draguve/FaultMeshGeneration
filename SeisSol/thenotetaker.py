from harness import *
import os
import re
import subprocess
from pprint import pprint
import numpy as np
import pandas as pd
import argparse
import sys

def get_all_sacct_info(job_id):
    try:
        # Construct the sacct command
        cmd = ["sacct", "-j", str(job_id), "--format=ALL", "-P"]
        # print(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        # print("\n--- STDOUT ---")
        # print(result.stdout)
        # print("\n--- STDERR ---")
        # print(result.stderr)
        if result.returncode != 0 or not result.stdout.strip():
            print("sacct failed or returned empty output.")
            return [{}]
        lines = result.stdout.strip().split("\n")
        # print(f"\nFound {len(lines)} lines.")
        headers = lines[0].split("|")
        # print(f"\nHeaders: {headers}")
        records = []
        for line in lines[1:]:
            # print(f"\nProcessing line: {line}")
            values = line.split("|")
            record = dict(zip(headers, values))
            records.append(record)
        # print(f"\nMatched {len(records)} entries.")
        return records
    except Exception as e:
        print(f"Error: {e}")
        return [{}]

def extract_key_value_pairs(file_path):
    # Define regex patterns for key-value pairs with colon or equals sign
    # The value must not contain line breaks (excluding spaces/tabs)
    patterns = [
        r'(\w+)\s*:\s*([\S]+)\n*$',  # Matches key: value (colon-based) where value doesn't have a newline
        r'(\w+)\s*=\s*([^\n;]+)'   # Matches key = value (equals-based) where value doesn't have a newline
    ]

    key_value_pairs = {}

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        content = file.read()

        # Search for both patterns
        for pattern in patterns:
            matches = re.findall(pattern, content,flags=re.MULTILINE)
            for match in matches:
                key, value = match
                key_value_pairs[key] = value.strip()  # Clean up any extra spaces in values
    return key_value_pairs

def extract_coordinates(expression):
    # Regular expression to match the x, y, z coordinates
    pattern = r'([xyz])([+-]?\d*\.\d+|\d+)'  # Match x, y, z followed by the number (with optional signs)
    # Find all matches for x, y, z with their corresponding values
    matches = re.findall(pattern, expression)
    # Create a dictionary to store the coordinates
    coordinates = {'x': None, 'y': None, 'z': None}
    # Assign the found values to the corresponding coordinate
    for axis, value in matches:
        coordinates[axis] = -float(value)
    # Return the extracted coordinates as a tuple
    return (coordinates['x'], coordinates['y'], coordinates['z'])


def get_job_id():
    # List files in the current directory
    files = os.listdir('.')
    # Regex pattern to match SLURM output filenames (e.g., 38483002.XY20-single-1m-ud0.3.out)
    pattern = re.compile(r'^(\d+)\..*\.out$')
    # Search for the job ID
    job_id = ""
    for f in files:
        match = pattern.match(f)
        if match:
            job_id = match.group(1)
            break  # remove this if you want to find all job IDs
    return job_id


def write_row(row):
    with open("notes.tsv", "w") as f:
        f.write("\t".join(map(str, row)) + "\n")


def get_max_values(params):
    output_prefix = params.get('OutputFile')
    energy_file = f"{output_prefix}-energy.csv"
    if not os.path.exists(energy_file):
        return "",""
    df = pd.read_csv(energy_file)
    df = df.pivot_table(index="time", columns="variable", values="measurement")
    df["M_0"] = df["seismic_moment"] - df["total_frictional_work"] - df["elastic_energy"] - df["elastic_kinetic_energy"]
    df["M_0"] = df["M_0"].where(df["M_0"] > 0)  # mask non-positive values
    df["M_w"] = (np.log10(df["M_0"]) - 9.1) / 1.5
    return f"{np.max(df['M_0']):.10e}",f"{np.max(df['M_w'])}"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--parameter',default="parameters.par")
    args = parser.parse_args()
    par_file = args.parameter

    if not os.path.exists(par_file):
        print("Cannot find parameter file",file=sys.stderr)
        exit(1)

    row = []
    params = parse_parameters(par_file)
    output_prefix = params.get('OutputFile')
    row.append(os.getcwd())

    rough_values = extract_key_value_pairs(params.get('ModelFileName'))
    print(rough_values)
    print(params)

    material_values = extract_key_value_pairs(params.get('MaterialFileName'))
    print(material_values)

    job_id_str = get_job_id()
    all_slurm_records = {}
    if job_id_str != "":
        job_id = int(get_job_id())
        all_slurm_records = get_all_sacct_info(job_id)
        if len(all_slurm_records) == 0:
            all_slurm_records = {}
        else:
            all_slurm_records = all_slurm_records[0]
    pprint(all_slurm_records)

    m_0,m_w = get_max_values(params)

    coordinates = extract_coordinates(rough_values.get("r"))
    row.append(rough_values.get("mu_d",""))
    row.append(rough_values.get("mu_s",""))
    row.append(rough_values.get("d_c",""))
    row.append(rough_values.get("s_xx",""))
    row.append(rough_values.get("s_yy",""))
    row.append(rough_values.get("s_zz",""))
    row.append(rough_values.get("s_xy",""))
    row.append(rough_values.get("r_crit",""))
    row.append(rough_values.get("Vs",""))
    row.append(material_values.get("rho",""))
    row.append(material_values.get("mu",""))
    row.append(material_values.get("lambda",""))
    row.append("") # For notes
    row.append(coordinates)
    row.append("single") # assume Single for now
    row.append(all_slurm_records.get("Elapsed",""))
    row.append(all_slurm_records.get("AllocNodes",""))
    row.append("CPU")  # assume CPU
    row.append(params.get("EndTime",""))
    row.append(m_0) # Max M_0
    row.append(m_w) # Max M_w
    row.append("") # Link
    row.append(job_id_str) # JobID
    row.append(all_slurm_records.get("Planned",""))
    row.append(all_slurm_records.get("Submit",""))
    row.append(all_slurm_records.get("NodeList",""))
    row.append(all_slurm_records.get("CPUTimeRAW",""))
    print(row)
    write_row(row)

if __name__ == '__main__':
    main()
