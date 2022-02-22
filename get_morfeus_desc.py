import numpy as np

from morfeus.sasa import SASA
from morfeus.dispersion import Dispersion
from morfeus.io import read_xyz
import os
import pandas as pd


def get_sasa(folder_path, filename):
    try:
        elements, coordinates = read_xyz(f'{os.path.join(folder_path, filename)}')

        sasa = SASA(elements, coordinates)
        sasa_array = np.array([sasa.atom_areas[i] for i in range(1, len(sasa.atom_areas)+1)])
    except Exception:
        return None

    return sasa_array


def get_pint(folder_path, filename):
    try:
        elements, coordinates = read_xyz(f'{os.path.join(folder_path, filename)}')

        disp = Dispersion(elements, coordinates)
        pint_array = np.array([disp.atom_p_int[i] for i in range(1, len(disp.atom_p_int)+1)])
    except Exception:
        return None 

    return pint_array


if __name__ == '__main__':
    df = pd.read_csv('reactants_cycloadd.csv')
    df['sasa'] = df['id'].apply(lambda x: get_sasa('xyz_files', f'geometry_{x}.xyz'))
    df['pint'] = df['id'].apply(lambda x: get_pint('xyz_files', f'geometry_{x}.xyz'))

    print(df.head())

    df.to_csv('morfeus_desc.csv')
    df.to_pickle('morfeus_desc.pkl') 