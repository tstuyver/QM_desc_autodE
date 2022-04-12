from argparse import ArgumentParser, Namespace
import os
import shutil
import autode as ade
import subprocess
import pandas as pd

from lib import create_logger
from lib import Molecule
from lib import get_opt_molecules
from lib import dft_scf

G16_PATH = "/home/gridsan/tstuyver/RMG_shared/Software/gaussian/g16"

parser = ArgumentParser()
parser.add_argument(
    "--ismiles", type=str, required=False, help="input smiles included in a .csv file"
)
parser.add_argument(
    "--solvent-name",
    type=str,
    default=None,
    help="whether or not to determine the descriptors in solvent",
)

# xtb optimization
parser.add_argument(
    "--xtb-folder", type=str, default="XTB_opt", help="folder for XTB optimization"
)
parser.add_argument(
    "--xtb-n-procs", type=int, default=20, help="number of process for optimization"
)

# DFT calculation
parser.add_argument(
    "--DFT-folder", type=str, default="DFT", help="folder for DFT calculation"
)
parser.add_argument(
    "--DFT-theory",
    type=str,
    default="b3lyp/def2svp",
    help="level of theory for the DFT calculation",
)
parser.add_argument(
    "--DFT-n-procs", type=int, default=20, help="number of process for DFT calculations"
)

args = parser.parse_args()

ade.Config.n_cores = args.xtb_n_procs
ade.Config.hmethod_conformers = False

name = os.path.splitext(args.ismiles)[0]
output_name = f"{name}_qm_descriptors"
logger = create_logger(name)

df = pd.read_csv(args.ismiles, index_col=0)

pwd = os.getcwd()

logger.info("starting structure generation")

if not os.path.isdir(args.xtb_folder):
    os.mkdir(args.xtb_folder)
os.chdir(args.xtb_folder)

molecule_list = [
    Molecule(x[0], x[1], args.solvent_name) for x in df[["id", "smiles"]].values
]

opt_molecules = get_opt_molecules(args, molecule_list, logger)

os.chdir(pwd)

# G16 DFT calculation
logger.info("starting Gaussian single-point calculations for the optimized geometries")
if not os.path.isdir(args.DFT_folder):
    os.mkdir(args.DFT_folder)

qm_descriptors_list = []
failed_molecules = []
for molecule in opt_molecules:
    try:
        shutil.copyfile(
            os.path.join(args.xtb_folder, molecule.xyz_file),
            os.path.join(args.DFT_folder, molecule.xyz_file),
        )
        qm_descriptors = dft_scf(
            args.DFT_folder,
            molecule,
            G16_PATH,
            args.DFT_theory,
            args.DFT_n_procs,
            logger,
            args.solvent_name,
        )
        qm_descriptors_list.append(qm_descriptors)
    except Exception as e:
        logger.error(
            f"Gaussian single-point calculations for {os.path.splitext(molecule.xyz_file)[0]} failed: {e}"
        )
        failed_molecules.append(molecule)
        try:
            for folder in [
                os.path.join(args.DFT_folder, "neutral/"),
                os.path.join(args.DFT_folder, "minus1/"),
                os.path.join(args.DFT_folder, "plus1/"),
                os.path.join(args.DFT_folder, "multiplicity/"),
            ]:
                for fname in os.listdir(folder):
                    if "core" in fname:
                        os.remove(os.path.join(folder, fname))
        except Exception:
            continue

df_qm_descriptors = pd.DataFrame(qm_descriptors_list)
df_qm_descriptors.to_csv(f"{output_name}.csv")
df_qm_descriptors.to_pickle(f"{output_name}.pkl")

#________________________________
# This might not be necessary

for molecule in failed_molecules:
    try:
        shutil.copyfile(
            os.path.join(args.xtb_folder, molecule.xyz_file),
            os.path.join(args.DFT_folder, molecule.xyz_file),
        )
        qm_descriptors = dft_scf(
            args.DFT_folder,
            molecule,
            G16_PATH,
            args.DFT_theory,
            args.DFT_n_procs,
            logger,
            args.solvent_name,
        )
        qm_descriptors_list.append(qm_descriptors)
    except Exception as e:
        logger.error(
            f"Gaussian single-point calculations for {os.path.splitext(molecule.xyz_file)[0]} failed: {e}"
        )
        try:
           for folder in [
                os.path.join(args.DFT_folder, "neutral/"),
                os.path.join(args.DFT_folder, "minus1/"),
                os.path.join(args.DFT_folder, "plus1/"),
                os.path.join(args.DFT_folder, "multiplicity/"),
            ]:
                for fname in os.listdir(folder):
                    if "core" in fname:
                        os.remove(os.path.join(folder, fname))
        except Exception:
            continue 

df_qm_descriptors = pd.DataFrame(qm_descriptors_list)
df_qm_descriptors.to_csv(f"{output_name}.csv")
df_qm_descriptors.to_pickle(f"{output_name}.pkl")

#________________________________

subprocess.run(
    [
        "tar",
        "cvzf",
        os.path.join(pwd, f"xtb_{name}.tar.gz"),
        os.path.join(pwd, args.xtb_folder),
    ]
)
subprocess.run(["rm", "-r", os.path.join(pwd, args.xtb_folder)])

subprocess.run(
    [
        "tar",
        "cvzf",
        os.path.join(pwd, f"dft_{name}.tar.gz"),
        os.path.join(pwd, args.DFT_folder),
    ]
)
subprocess.run(["rm", "-r", os.path.join(pwd, args.DFT_folder)])
