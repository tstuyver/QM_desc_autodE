from argparse import ArgumentParser, Namespace
import os
import shutil
import autode as ade
import subprocess
import pandas as pd

from lib import create_logger
from lib import Molecule
from lib import get_opt_molecules
from lib import dft_scf, extract_descriptors
from lib import read_log

G16_PATH = "XXX"

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
parser.add_argument(
    "--resume-fail",
    action="store_true",
    help="whether to just resume in case of previous failure",
)

# xtb optimization
parser.add_argument(
    "--xtb-folder", type=str, default="XTB_opt", help="folder for XTB optimization"
)
parser.add_argument(
    "--xtb-n-procs", type=int, default=48, help="number of process for optimization"
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
    "--DFT-n-procs", type=int, default=48, help="number of process for DFT calculations"
)


def get_qm_descriptors(args, opt_molecules):
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
                os.remove(os.path.join(args.DFT_folder, molecule.xyz_file))
                continue

    return qm_descriptors_list, failed_molecules


if __name__ == "__main__":
    args = parser.parse_args()

    # set autodE variables
    ade.Config.n_cores = args.xtb_n_procs
    ade.Config.hmethod_conformers = False

    name = os.path.splitext(args.ismiles)[0]
    output_name = f"{name}_qm_descriptors"
    logger = create_logger(name)

    df = pd.read_csv(args.ismiles, index_col=0)

    pwd = os.getcwd()

    # start from scratch
    if not args.resume_fail:
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
        logger.info(
            "starting Gaussian single-point calculations for the optimized geometries"
        )
        if not os.path.isdir(args.DFT_folder):
            os.mkdir(args.DFT_folder)

        qm_descriptors_list, failed_molecules = get_qm_descriptors(args, opt_molecules)

    # resume after time out or failure
    if args.resume_fail:
        logger.info("resuming descriptor computation")
        # print(os.listdir(args.xtb_folder))
        opt_molecules = [
            os.path.splitext(os.path.basename(file_name))[0]
            for file_name in os.listdir(args.xtb_folder)
            if file_name.endswith(".xyz")
        ]

        # print(opt_molecules)
        opt_molecules.sort()
        qm_descriptors_list = []
        ids_to_do = []
        os.chdir(args.DFT_folder)

        for molecule in opt_molecules:
            qm_descriptors = extract_descriptors(molecule, df)

            if qm_descriptors is not None:
                qm_descriptors_list.append(qm_descriptors)
            else:
                ids_to_do.append(molecule)

        os.chdir(pwd)
        molecule_list = []

        for id in ids_to_do:
            molecule_list.append(
                Molecule(
                    df.iloc[int(id.split("_")[-1])]["id"],
                    df.iloc[int(id.split("_")[-1])]["smiles"],
                    args.solvent_name,
                )
            )

        for molecule in molecule_list:
            molecule.xyz_file = f"geometry_{molecule.id}.xyz"

        qm_descriptors_list_tmp, failed_molecules = get_qm_descriptors(
            args, molecule_list
        )
        for qm_descriptors in qm_descriptors_list_tmp:
            qm_descriptors_list.append(qm_descriptors)

    df_qm_descriptors = pd.DataFrame(qm_descriptors_list)
    df_qm_descriptors.to_csv(f"{output_name}.csv")
    df_qm_descriptors.to_pickle(f"{output_name}.pkl")

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
