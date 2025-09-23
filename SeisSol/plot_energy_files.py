import glob
import os
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import typer
from typing_extensions import Annotated, Optional


def get_units(variable_name):
    units = {
        "gravitational_energy": "J",
        "acoustic_energy": "J",
        "acoustic_kinetic_energy": "J",
        "elastic_energy": "J",
        "elastic_kinetic_energy": "J",
        "momentumX": "kg·m/s",
        "momentumY": "kg·m/s",
        "momentumZ": "kg·m/s",
        "total_frictional_work": "J",
        "static_frictional_work": "J",
        "seismic_moment": "N·m",
        "potency": "m³",
        "plastic_moment": "J",
        "seismic_moment_rate": "N·m/s",
        "M_0": "N·m",
        "M_w": "",
    }
    return units.get(variable_name, "")


def convert_to_df(filepath: Path) -> pd.DataFrame:
    if not os.path.exists(filepath):
        print(f"Could not find energy {filepath}!")
    df = pd.read_csv(filepath)
    df = df.pivot_table(index="time", columns="variable", values="measurement")
    df["seismic_moment_rate"] = np.gradient(df["seismic_moment"], df.index[1])
    df["M_0"] = (
        df["seismic_moment"]
        - df["total_frictional_work"]
        - df["elastic_energy"]
        - df["elastic_kinetic_energy"]
    )
    df["M_0"] = df["M_0"].where(df["M_0"] > 0)  # mask non-positive values
    df["M_w"] = (np.log10(df["M_0"]) - 9.1) / 1.5
    return df


def plot_and_save_column(
    column_name: str,
    file_data_frames: list[tuple[str, pd.DataFrame]],
    output_path: Path,
):
    print(f"Plotting {column_name}")
    plt.xlabel("Time (s)")
    plt.ylabel(f"{column_name} in {get_units(column_name)}")
    for filename, df in file_data_frames:
        if column_name in df.columns:
            time = df.index.to_numpy()
            y_axis = df[column_name].to_numpy()
            plt.plot(time, y_axis, label=filename)
    plt.legend()
    plt.savefig(output_path)
    plt.clf()


def main(
    energy_directory: Annotated[
        Path, typer.Argument(help="Directory containing all the energy files")
    ],
    output: Annotated[
        Optional[Path],
        typer.Option(
            help="Output path for the plots, will make direcory if it does not exist, by default will add plots folder to energy direcory"
        ),
    ] = None,
):
    if not os.path.isdir(energy_directory):
        print(f"Could not find {energy_directory}!")
        exit(1)

    if output is None:
        output = energy_directory / "plots/"

    output.mkdir(parents=True, exist_ok=True)

    csv_files = glob.glob(os.path.join(energy_directory, "*.csv"))
    file_data_frames = [
        (Path(file).stem, convert_to_df(Path(file))) for file in csv_files
    ]
    all_columns = set.union(*[set(df.columns) for _, df in file_data_frames])
    for column in all_columns:
        output_path = output / f"{column}.svg"
        plot_and_save_column(column, file_data_frames, output_path)


if __name__ == "__main__":
    typer.run(main)
