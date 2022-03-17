from rdkit import Chem
from rdkit.Chem import Descriptors
import os
import subprocess
from .file_parser import xyz2com
from .grab_QM_descriptors import read_log

solvent_dict = {"water": "Water", "dmf": "n,n-DiMethylFormamide"}


def dft_scf(folder, molecule, g16_path, level_of_theory, n_procs, logger, solvent_name):
    basename = os.path.basename(molecule.xyz_file)

    parent_folder = os.getcwd()
    os.chdir(folder)

    try:
        file_name = os.path.splitext(basename)[0]
        xc = level_of_theory.split("/")[0]
        ao = level_of_theory.split("/")[1]

        mol = Chem.MolFromSmiles(molecule.smiles, sanitize=False)
        mol_no_H = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(mol_no_H)
        mol_charge = Chem.GetFormalCharge(mol_no_H)
        mol_radical_electrons = Descriptors.NumRadicalElectrons(mol_no_H)

        with open(f"{molecule.xyz_file}", "r") as f:
            xyz = ""
            for line in f.readlines():
                xyz += line

        pwd = os.getcwd()

        g16_command = os.path.join(g16_path, "g16")
        QM_descriptors = {}
        for jobtype in ["neutral", "plus1", "minus1", "multiplicity"]:
            if not os.path.isdir(jobtype):
                os.mkdir(jobtype)

            if jobtype == "neutral":
                charge = mol_charge
                if mol_radical_electrons == 0:
                    mult = 1
                elif mol_radical_electrons == 1:
                    mult = 2
                if solvent_name:
                    head = (
                        f"%chk={file_name}.chk\n%nprocshared={n_procs}\n# {xc}/{ao} nmr=GIAO scf=(maxcycle=512, xqc) "
                        f"pop=(full,mbs,hirshfeld,nboread) scrf=(SMD, Solvent={solvent_dict[solvent_name]})\n"
                    )
                else:
                    head = (
                        f"%chk={file_name}.chk\n%nprocshared={n_procs}\n# {xc}/{ao} nmr=GIAO scf=(maxcycle=512, xqc) "
                        "pop=(full,mbs,hirshfeld,nboread)\n"
                    )
            elif jobtype == "plus1":
                charge = mol_charge + 1
                if mol_radical_electrons == 0:
                    mult = 2
                elif mol_radical_electrons == 1:
                    mult = 1
                head = (
                    "%chk={}.chk\n%nprocshared={}\n# {}/{} scf=(maxcycle=512, xqc) "
                    "pop=(full,mbs,hirshfeld,nboread)\n".format(
                        file_name, n_procs, xc, ao
                    )
                )
            elif jobtype == "minus1":
                charge = mol_charge - 1
                if mol_radical_electrons == 0:
                    mult = 2
                elif mol_radical_electrons == 1:
                    mult = 1
                head = (
                    "%chk={}.chk\n%nprocshared={}\n# {}/{} scf=(maxcycle=512, xqc) "
                    "pop=(full,mbs,hirshfeld,nboread)\n".format(
                        file_name, n_procs, xc, ao
                    )
                )
            elif jobtype == "multiplicity":
                charge = mol_charge
                if mol_radical_electrons == 0:
                    mult = 3
                elif mol_radical_electrons == 1:
                    mult = 2
                head = (
                    "%chk={}.chk\n%nprocshared={}\n# {}/{} scf=(maxcycle=512, xqc) "
                    "pop=(full,mbs,hirshfeld,nboread)\n".format(
                        file_name, n_procs, xc, ao
                    )
                )

            os.chdir(jobtype)
            comfile = file_name + ".gjf"
            xyz2com(
                xyz,
                head=head,
                comfile=comfile,
                charge=charge,
                mult=mult,
                footer="$NBO BNDIDX $END\n",
            )

            logfile = file_name + ".log"
            outfile = file_name + ".out"
            with open(outfile, "w") as out:
                subprocess.run(
                    "{} < {} >> {}".format(g16_command, comfile, logfile),
                    shell=True,
                    stdout=out,
                    stderr=out,
                )
                QM_descriptors[jobtype] = read_log(logfile, jobtype, smiles)
            os.chdir(pwd)

        QM_descriptor_return = QM_descriptors["neutral"]

        # charges and fukui indices
        for charge in ["mulliken_charge", "hirshfeld_charges", "NPA_Charge"]:
            QM_descriptor_return["{}_plus1".format(charge)] = QM_descriptors["plus1"][
                charge
            ]
            QM_descriptor_return["{}_minus1".format(charge)] = QM_descriptors["minus1"][
                charge
            ]

            QM_descriptor_return["{}_fukui_elec".format(charge)] = (
                QM_descriptors["neutral"][charge] - QM_descriptors["minus1"][charge]
            )
            QM_descriptor_return["{}_fukui_neu".format(charge)] = (
                QM_descriptors["plus1"][charge] - QM_descriptors["neutral"][charge]
            )

            # spin density
        for spin in [
            "mulliken_spin_density",
            "hirshfeld_spin_density",
            "NPA_spin_density",
        ]:
            QM_descriptor_return["{}_plus1".format(spin)] = QM_descriptors["plus1"][
                spin
            ]
            QM_descriptor_return["{}_minus1".format(spin)] = QM_descriptors["minus1"][
                spin
            ]
            QM_descriptor_return["{}_multiplicity".format(spin)] = QM_descriptors[
                "multiplicity"
            ][spin]

        # SCF
        QM_descriptor_return["SCF_plus1"] = QM_descriptors["plus1"]["SCF"]
        QM_descriptor_return["SCF_minus1"] = QM_descriptors["minus1"]["SCF"]
        QM_descriptor_return["SCF_multiplicity"] = QM_descriptors["multiplicity"]["SCF"]

        os.remove(molecule.xyz_file)
        logger.info(
            "Gaussian single-point calculations for {} completed".format(file_name)
        )
    finally:
        os.chdir(parent_folder)

    return QM_descriptor_return
