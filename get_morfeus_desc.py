import numpy as np

from morfeus.sasa import SASA
from morfeus.dispersion import Dispersion
from morfeus.io import read_xyz
import os
import pandas as pd
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "--csv-file",
    type=str,
    required=True,
    help="input .csv file containing list of molecules",
)
parser.add_argument(
    "--xyz-folder",
    type=str,
    required=True,
    help="folder containing .xyz files for each molecule",
)


def get_sasa(folder_path, filename):
    try:
        elements, coordinates = read_xyz(f"{os.path.join(folder_path, filename)}")

        sasa = SASA(elements, coordinates)
        sasa_array = np.array(
            [sasa.atom_areas[i] for i in range(1, len(sasa.atom_areas) + 1)]
        )
    except Exception:
        return None

    return sasa_array


def get_pint(folder_path, filename):
    try:
        elements, coordinates = read_xyz(f"{os.path.join(folder_path, filename)}")

        disp = Dispersion(elements, coordinates)
        pint_array = np.array(
            [disp.atom_p_int[i] for i in range(1, len(disp.atom_p_int) + 1)]
        )
    except Exception:
        return None

    return pint_array


if __name__ == "__main__":
    args = parser.parse_args()
    path = "/".join(args.csv_file.split("/")[:-1])
    df = pd.read_csv(args.csv_file)
    df["sasa"] = df["id"].apply(
        lambda x: get_sasa(args.xyz_folder, f"geometry_{x}.xyz")
    )
    df["pint"] = df["id"].apply(
        lambda x: get_pint(args.xyz_folder, f"geometry_{x}.xyz")
    )

    print(df.head())

    df.to_csv(f"{path}/morfeus_desc.csv")
    df.to_pickle(f"{path}/morfeus_desc.pkl")
