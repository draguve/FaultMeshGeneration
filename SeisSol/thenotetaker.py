from harness import *
import os
import re
import subprocess
from datetime import datetime

def get_slurm_job_usage(job_id):
    try:
        # Request specific fields from sacct
        result = subprocess.run(
            [
                "sacct",
                "-j", str(job_id),
                "--format=JobIDRaw,Start,End,Elapsed,NCPUS,NNodes",
                "--noheader",
                "--parsable2"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"sacct error: {result.stderr.strip()}")
        # Find the line corresponding to the base job ID (not .batch or .extern)
        lines = result.stdout.strip().split("\n")
        for line in lines:
            parts = line.strip().split("|")
            if parts[0] == str(job_id):
                start, end, elapsed, ncpus, nnodes = parts[1], parts[2], parts[3], parts[4], parts[5]
                return elapsed, int(ncpus), int(nnodes)
        return None, None, None
    except Exception as e:
        print(f"Error: {e}")
        return None, None, None

def extract_key_value_pairs(file_path):
    # Define regex patterns for key-value pairs with colon or equals sign
    # The value must not contain line breaks (excluding spaces/tabs)
    patterns = [
        r'(\w+)\s*:\s*([\S]+)\n',  # Matches key: value (colon-based) where value doesn't have a newline
        r'(\w+)\s*=\s*([^\n;]+)'   # Matches key = value (equals-based) where value doesn't have a newline
    ]

    key_value_pairs = {}

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        content = file.read()

        # Search for both patterns
        for pattern in patterns:
            matches = re.findall(pattern, content)
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
    if not os.path.exists("parameters.par"):
        print("Cannot find parameters.par file",file=sys.stderr)
        exit(1)

    row = []
    params = parse_parameters("parameters.par")
    output_prefix = params.get('OutputFile')
    row.append(os.getcwd())

    rough_values = extract_key_value_pairs(params.get('ModelFileName'))
    print(rough_values)
    print(params)

    job_id = int(get_job_id())
    elapsed, cores, nodes = get_slurm_job_usage(job_id)

    m_0,m_w = get_max_values(params)

    coordinates = extract_coordinates(rough_values.get("r"))
    row.append(rough_values.get("mu_d"))
    row.append(rough_values.get("mu_s"))
    row.append(rough_values.get("d_c"))
    row.append(rough_values.get("s_xx"))
    row.append(rough_values.get("s_yy"))
    row.append(rough_values.get("s_zz"))
    row.append(rough_values.get("s_xy"))
    row.append(rough_values.get("r_crit"))
    row.append(rough_values.get("Vs"))
    row.append("") # For notes
    row.append(coordinates)
    row.append("single") # assume Single for now
    row.append(elapsed) #Exec time
    row.append(nodes) # assume 4 nodes for now
    row.append("CPU")  # assume CPU
    row.append(params.get("EndTime"))
    row.append(m_0) # Max M_0
    row.append(m_w) # Max M_w
    row.append("") # Link
    row.append(get_job_id()) # JobID
    write_row(row)




if __name__ == '__main__':
    main()