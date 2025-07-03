import os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import argparse
from pathlib import Path
import sys

def parse_parameters(file_path):
    output = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Remove comments
            line = line.split('!')[0].strip()
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip().strip("'\"")  # Remove quotes if any
                output[key] = value
    return output

def plot_and_save(df,what_to_plot,plots_path):
    df.plot(y=what_to_plot,use_index=True)
    plt.savefig(f"{plots_path}{what_to_plot}.png", dpi=300, bbox_inches='tight')

def before(params):
    output_prefix = os.path.abspath(params.get('OutputFile'))
    output_dir = os.path.dirname(output_prefix)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(output_dir)

def after(params):
    if int(params.get("EnergyOutput")) != 1:
        print("Energy Output is disabled")
        exit(0)
    output_prefix = params.get('OutputFile')
    plots_path = f"{output_prefix}-plots/"
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)
    energy_file = f"{output_prefix}-energy.csv"
    if not os.path.exists(energy_file):
        print("Could not find energy file!")
        exit(1)
    print("Generating plots...")
    df = pd.read_csv(energy_file)
    df = df.pivot_table(index="time", columns="variable", values="measurement")
    df["seismic_moment_rate"] = np.gradient(df["seismic_moment"], df.index[1]) #df index gets the time difference between 2 points
    plot_and_save(df,"seismic_moment_rate",plots_path)
    df["M_0"] = df["seismic_moment"] - df["total_frictional_work"] - df["elastic_energy"] - df["elastic_kinetic_energy"]
    df["M_0"] = df["M_0"].where(df["M_0"] > 0)  # mask non-positive values
    df["M_w"] = (np.log10(df["M_0"]) - 9.1) / 1.5

    plot_and_save(df,"M_0",plots_path)
    plot_and_save(df,"M_w",plots_path)

    plot_and_save(df,"acoustic_energy",plots_path)
    plot_and_save(df,"acoustic_kinetic_energy",plots_path)
    plot_and_save(df,"elastic_energy",plots_path)
    plot_and_save(df,"elastic_kinetic_energy",plots_path)
    plot_and_save(df,"gravitational_energy",plots_path)
    plot_and_save(df,"momentumX",plots_path)
    plot_and_save(df,"momentumY",plots_path)
    plot_and_save(df,"momentumZ",plots_path)
    plot_and_save(df,"plastic_moment",plots_path)
    plot_and_save(df,"potency",plots_path)
    plot_and_save(df,"seismic_moment",plots_path)
    plot_and_save(df,"static_frictional_work",plots_path)
    plot_and_save(df,"total_frictional_work",plots_path)

    print(f"Max M_0 : {np.max(df['M_0']):.10e}")
    print(f"Max M_w : {np.max(df['M_w'])}")

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'mode',
        choices=['before', 'after'],
        help="Mode must be 'before' or 'after'"
    )
    parser.add_argument('-p', '--parameter',default="parameters.par")
    args = parser.parse_args()

    par_file = args.parameter

    if not os.path.exists(par_file):
        print("Cannot find parameters.par file",file=sys.stderr)
        exit(1)
    params = parse_parameters(par_file)

    if args.mode == "before":
        before(params)
    elif args.mode == "after":
        after(params)

if __name__ == '__main__':
    main()
