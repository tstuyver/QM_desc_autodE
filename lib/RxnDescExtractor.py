import pandas as pd
import rdkit.Chem as Chem
import numpy as np

hartree = 627.5094740631


def strip_atom_map_num(smiles):
    """Strip atom map numbers from SMILES.

    Args:
        smiles (str): The initial atom-mapped SMILES string

    Returns:
        str: SMILES string with map numbers removed
    """
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    stripped_smiles = Chem.MolToSmiles(mol)

    return stripped_smiles


def turn_elements_into_arrays(df, column_name):
    """Formats the elements in a column of a dataframe as array
    (needed due to incompatibility between pandas versions for morfeus-descriptors).

    Args:
        df (pd.DataFrame): unformatted dataframe
        column_name (str): name of the column

    Returns:
        pd. DataFrame: formatted dataframe
    """
    df[f"{column_name}"] = df[f"{column_name}"].apply(lambda x: np.array([float(x)]))
    return df


class RxnDescExtractor:
    """Class for reaction descriptor extraction.

    Attributes:
        desc_file (str): path to the descriptor .pkl file
        reactions_file (str): path to the dataframe containing all the reaction SMILES
        output_name (str): name of the output_file to be created
    """

    def __init__(self, desc_file, reactions_file, output_name):
        self.df_desc = pd.read_pickle(desc_file)
        self.df_reactions = pd.read_csv(reactions_file)
        self.output_name = output_name

        self.mol_desc = self.extract_descriptor_values()

        self.reaction_desc = self.determine_reaction_descs()

    def return_valid_df_reactions(self):
        """Drop invalid reaction SMILES"""
        # df_return = self.df_reactions.dropna()[['rxn_id', 'rxn_smiles', 'solvent', 'temp', 'G_act', 'G_r']]
        df_return = self.df_reactions.dropna()[["rxn_id", "rxn_smiles"]]
        return df_return

    def determine_reaction_descs(self):
        """Determine the reaction descriptors. Several versions of promotion gaps are computed:
        1) gaps based on HOMO/LUMO energies (orb)
        2) gaps based on SCF energies of cation/anion
        3) gaps based on SCF energies of cation/anion, and corrected according to the scheme by De Proft and Tozer (https://doi.org/10.1039/B605302P)

        Additionally, also orbital energies are computed.

        """
        self.df_reactions["G_orb"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_orbitals(x)
        )
        self.df_reactions["G"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_scf(x)
        )
        self.df_reactions["G_alt1"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_alt_scf(x, 1, corr=False)
        )
        self.df_reactions["G_alt2"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_alt_scf(x, 2, corr=False)
        )

        self.df_reactions["G_alt1_orb"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_alt_orbitals(x, 1)
        )
        self.df_reactions["G_alt1_corr"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_alt_scf(x, 1)
        )
        self.df_reactions["G_alt2_orb"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_alt_orbitals(x, 2)
        )
        self.df_reactions["G_alt2_corr"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_promotion_gap_alt_scf(x, 2)
        )

        self.df_reactions["orbital_energies"] = self.df_reactions["rxn_smiles"].apply(
            lambda x: self.get_orbital_energies(x)
        )
        self.df_reactions["homo_1"] = self.df_reactions["orbital_energies"].apply(
            lambda x: x[0]
        )
        self.df_reactions["lumo_1"] = self.df_reactions["orbital_energies"].apply(
            lambda x: x[1]
        )
        self.df_reactions["homo_2"] = self.df_reactions["orbital_energies"].apply(
            lambda x: x[2]
        )
        self.df_reactions["lumo_2"] = self.df_reactions["orbital_energies"].apply(
            lambda x: x[3]
        )

    def extract_descriptor_values(self):
        """Extract descriptor values and store them in separate dicts."""
        mol_desc = {
            "reactants_singlet": dict(zip(self.df_desc["smiles"], self.df_desc["SCF"])),
            "reactants_homo": dict(zip(self.df_desc["smiles"], self.df_desc["homo"])),
            "reactants_lumo": dict(zip(self.df_desc["smiles"], self.df_desc["lumo"])),
            "reactants_triplet": dict(
                zip(self.df_desc["smiles"], self.df_desc["SCF_multiplicity"])
            ),
            "reactants_plus": dict(
                zip(self.df_desc["smiles"], self.df_desc["SCF_plus1"])
            ),
            "reactants_minus": dict(
                zip(self.df_desc["smiles"], self.df_desc["SCF_minus1"])
            ),
        }

        return mol_desc

    def get_orbital_energies(self, rxn_smiles):
        """Get the orbital energies for an individual reaction SMILES."""
        reactants = rxn_smiles.split(">")[0].split(".")
        orbital_energies = []
        for reactant in reactants:
            try:
                stripped_reactant = strip_atom_map_num(reactant)
                orbital_energies.append(
                    self.mol_desc["reactants_homo"][stripped_reactant]
                )
                orbital_energies.append(
                    self.mol_desc["reactants_lumo"][stripped_reactant]
                )
            except Exception as e:
                print(f"Error for {rxn_smiles}:{e}")
                return [None, None, None, None]

        return orbital_energies

    def get_promotion_gap_orbitals(self, rxn_smiles):
        """Get the orbital-based promotion gap for an individual reaction SMILES."""
        reactants = rxn_smiles.split(">")[0].split(".")
        promotion_energy = 0
        for reactant in reactants:
            try:
                stripped_reactant = strip_atom_map_num(reactant)
                promotion_energy += (
                    self.mol_desc["reactants_lumo"][stripped_reactant]
                    - self.mol_desc["reactants_homo"][stripped_reactant]
                ) * hartree
            except Exception as e:
                print(f"Error for {rxn_smiles}: {e}")
                promotion_energy = None

        return promotion_energy

    def get_promotion_gap_scf(self, rxn_smiles):
        """Get the SCF-based promotion gap for an individual reaction SMILES."""
        reactants = rxn_smiles.split(">")[0].split(".")
        promotion_energy = 0
        for reactant in reactants:
            try:
                stripped_reactant = strip_atom_map_num(reactant)
                promotion_energy += (
                    self.mol_desc["reactants_triplet"][stripped_reactant]
                    - self.mol_desc["reactants_singlet"][stripped_reactant]
                ) * hartree
            except Exception as e:
                print(f"Error for {rxn_smiles}: {e}")
                promotion_energy = None

        return promotion_energy

    def get_promotion_gap_alt_orbitals(self, rxn_smiles, num):
        """Get the orbital-based alternative promotion gap for an individual reaction SMILES."""
        reactants = rxn_smiles.split(">")[0].split(".")
        for reactant in reactants:
            reactant_mol = Chem.MolFromSmiles(reactant)
            try:
                if num == 1:
                    if any(
                        [
                            atom.GetFormalCharge() == 1
                            for atom in reactant_mol.GetAtoms()
                        ]
                    ):
                        dipole = strip_atom_map_num(reactant)
                        electron_affinity = (
                            -self.mol_desc["reactants_lumo"][dipole] * hartree
                        )
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        ionisation_potential = (
                            -self.mol_desc["reactants_homo"][dipolarophile] * hartree
                        )
                if num == 2:
                    if any(
                        [
                            atom.GetFormalCharge() == 1
                            for atom in reactant_mol.GetAtoms()
                        ]
                    ):
                        dipole = strip_atom_map_num(reactant)
                        ionisation_potential = (
                            -self.mol_desc["reactants_homo"][dipole] * hartree
                        )
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        electron_affinity = (
                            -self.mol_desc["reactants_lumo"][dipolarophile] * hartree
                        )
            except Exception as e:
                print(f"Error for {rxn_smiles}: {e}")
                return None

        return ionisation_potential - electron_affinity

    def get_promotion_gap_alt_scf(self, rxn_smiles, num, corr=True):
        """Get the SCF-based alternative promotion gap for an individual reaction SMILES."""
        reactants = rxn_smiles.split(">")[0].split(".")
        for reactant in reactants:
            reactant_mol = Chem.MolFromSmiles(reactant)
            try:
                if num == 1:
                    if any(
                        [
                            atom.GetFormalCharge() == 1
                            for atom in reactant_mol.GetAtoms()
                        ]
                    ):
                        dipole = strip_atom_map_num(reactant)
                        electron_affinity = (
                            self.mol_desc["reactants_singlet"][dipole]
                            - self.mol_desc["reactants_minus"][dipole]
                        ) * hartree
                        if electron_affinity < 0 and corr:
                            electron_affinity = (
                                -(
                                    self.mol_desc["reactants_lumo"][dipole]
                                    + self.mol_desc["reactants_homo"][dipole]
                                )
                                - (
                                    self.mol_desc["reactants_plus"][dipole]
                                    - self.mol_desc["reactants_singlet"][dipole]
                                )
                            ) * hartree
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        ionisation_potential = (
                            self.mol_desc["reactants_plus"][dipolarophile]
                            - self.mol_desc["reactants_singlet"][dipolarophile]
                        ) * hartree
                if num == 2:
                    if any(
                        [
                            atom.GetFormalCharge() == 1
                            for atom in reactant_mol.GetAtoms()
                        ]
                    ):
                        dipole = strip_atom_map_num(reactant)
                        ionisation_potential = (
                            self.mol_desc["reactants_plus"][dipole]
                            - self.mol_desc["reactants_singlet"][dipole]
                        ) * hartree
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        electron_affinity = (
                            self.mol_desc["reactants_singlet"][dipolarophile]
                            - self.mol_desc["reactants_minus"][dipolarophile]
                        ) * hartree
                        if electron_affinity < 0 and corr:
                            electron_affinity = (
                                -(
                                    self.mol_desc["reactants_lumo"][dipolarophile]
                                    + self.mol_desc["reactants_homo"][dipolarophile]
                                )
                                - (
                                    self.mol_desc["reactants_plus"][dipolarophile]
                                    - self.mol_desc["reactants_singlet"][dipolarophile]
                                )
                            ) * hartree
            except Exception as e:
                print(f"Error for {rxn_smiles}: {e}")
                return None

        return ionisation_potential - electron_affinity

    def output_reaction_descs_chemprop(self, descriptor_list=["G", "G_alt1", "G_alt2"]):
        """Output the reaction descriptors in chemprop format."""
        df_chemprop = self.df_reactions[[name for name in descriptor_list]]
        df_chemprop = df_chemprop.dropna()
        for name in descriptor_list:
            df_chemprop = turn_elements_into_arrays(df_chemprop, name)

        df_chemprop.to_pickle(
            f"reaction_desc_{self.output_name}_chemprop.pkl", protocol=4
        )
        df_chemprop.to_csv(f"reaction_desc_{self.output_name}_chemprop.csv")

    def output_reaction_descs_wln(self, descriptor_list=["G", "G_alt1", "G_alt2"]):
        """Output the reaction descriptors in wln format."""
        df_wln = self.df_reactions[["rxn_smiles"] + [name for name in descriptor_list]]
        df_wln = df_wln.rename(columns={"rxn_smiles": "smiles"})
        df_wln.dropna()
        df_wln = df_wln.drop_duplicates(subset=["smiles"])

        df_wln.to_pickle(f"reaction_desc_{self.output_name}_wln.pkl", protocol=4)
        df_wln.to_csv(f"reaction_desc_{self.output_name}_wln.csv")
