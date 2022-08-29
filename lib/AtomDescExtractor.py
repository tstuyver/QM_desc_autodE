import pandas as pd
import rdkit.Chem as Chem
import numpy as np

hartree = 627.5094740631


def strip_atom_map_num(smiles):
    """Strip atom map numbers from SMILES

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


def verify_smiles(numbered_rxn_smiles, compound_desc_list):
    """Check if all reactant SMILES in a reaction are present in the compound_des_list

    Args:
        numbered_rxn_smiles (str): input reaction SMILES
        compound_desc_list (List[str]): list of SMILES strings

    Returns:
        bool: True (all SMILES present in compound_desc_list) or False (when not the case)
    """
    for smiles in numbered_rxn_smiles.split(">")[0].split("."):
        if strip_atom_map_num(smiles) in compound_desc_list:
            continue
        else:
            return False

    return True


def reorder_entry(entry_wrong_order, ordering):
    """Reorder a descriptor array according to the atom-mapped SMILES.

    Args:
        entry_wrong_order (np.array): original descriptor array where atoms are not ordered in the right way
        ordering (List[int]): correct position for each of the indices in the original array

    Returns:
        np.array: the descriptor array in the right order
    """
    try:
        entry_correct_order = [0] * len(entry_wrong_order)
        for i in range(len(ordering)):
            entry_correct_order[ordering[i]] = entry_wrong_order[i]
        for j in range(len(ordering), len(entry_wrong_order)):
            entry_correct_order[j] = entry_wrong_order[j]

        return np.array(entry_correct_order)
    except IndexError:
        return None


def turn_elements_into_arrays(df):
    """Formats the elements in a column of a dataframe as array
    (needed due to incompatibility between pandas versions for morfeus-descriptors).

    Args:
        df (pd.DataFrame): unformatted dataframe
        column_name (str): name of the column

    Returns:
        pd. DataFrame: formatted dataframe
    """
    df = df.apply(lambda x: np.array(list(map(float, x))))

    return df


def match_smiles(smiles, numbered_smiles_list):
    """Match smiles with numbered_smiles_list."""
    for numbered_smiles in numbered_smiles_list:
        if smiles == strip_atom_map_num(numbered_smiles):
            return numbered_smiles


def get_bond_orders(smiles, bond_orders, ordering):
    """Construct and array of bond orders in the correct order based on the SMILES."""
    bond_order_final = []
    bond_order_dict = {}

    mol = Chem.MolFromSmiles(strip_atom_map_num(smiles))
    mol = Chem.AddHs(mol, addCoords=True)

    # construct bond_order_dict
    for bond in mol.GetBonds():
        if (
            bond.GetBeginAtom().GetSymbol() != "H"
            and bond.GetEndAtom().GetSymbol() != "H"
        ):
            i, j = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
            bond_order_dict[frozenset({i, j})] = bond_orders[i, j]

    mol_no_H = Chem.MolFromSmiles(smiles)

    # append the bond order values in the right order
    for bond in mol_no_H.GetBonds():
        if (
            bond.GetBeginAtom().GetSymbol() != "H"
            and bond.GetEndAtom().GetSymbol() != "H"
        ):
            i, j = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
            bond_order_final.append(
                bond_order_dict[frozenset({ordering.index(i), ordering.index(j)})]
            )

    return np.array(bond_order_final)


def get_bond_lengths(smiles, coord_list, ordering):
    """Construct an array of bond lengths in the correct order based on the SMILES"""
    bond_lengths = []
    bond_length_dict = {}
    distance_matrix = get_distance_matrix(coord_list)

    mol = Chem.MolFromSmiles(strip_atom_map_num(smiles))
    mol = Chem.AddHs(mol, addCoords=True)

    # construct bond_length_dict
    for bond in mol.GetBonds():
        if (
            bond.GetBeginAtom().GetSymbol() != "H"
            and bond.GetEndAtom().GetSymbol() != "H"
        ):
            i, j = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
            bond_length_dict[frozenset({i, j})] = distance_matrix[i, j]

    mol_no_H = Chem.MolFromSmiles(smiles)

    # append the bond length values in the right order
    for bond in mol_no_H.GetBonds():
        if (
            bond.GetBeginAtom().GetSymbol() != "H"
            and bond.GetEndAtom().GetSymbol() != "H"
        ):
            i, j = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
            bond_lengths.append(
                bond_length_dict[frozenset({ordering.index(i), ordering.index(j)})]
            )

    return np.array(bond_lengths)


def get_distance_matrix(coord_list):
    """Get distances from the coordinates."""
    distance_matrix = np.zeros((len(coord_list), len(coord_list)))
    for i in range(len(coord_list)):
        for j in range(len(coord_list)):
            distance_matrix[i][j] = np.sqrt(
                (coord_list[i][0] - coord_list[j][0]) ** 2
                + (coord_list[i][1] - coord_list[j][1]) ** 2
                + (coord_list[i][2] - coord_list[j][2]) ** 2
            )

    return distance_matrix


def include_morfeus_data(smiles, morfeus_dict):
    """Get morfeus descriptors."""
    try:
        return morfeus_dict[smiles]
    except KeyError:
        return None


class AtomDescExtractor:
    """Class for atom descriptor extraction.

    Attributes:
        desc_file (str): path to the descriptor .pkl file
        reactions_file (str): path to the dataframe containing all the reaction SMILES
        output_name (str): name of the output_file to be created
        morfeus_file (str): path to the morfeus descriptor file (default: None)
        partitioning_scheme (str): which partitioning scheme has been used (default: Hirshfeld)
    """

    def __init__(
        self,
        desc_file,
        reactions_file,
        output_name,
        morfeus_file=None,
        partitioning_scheme="hirshfeld",
    ):
        self.df_desc = pd.read_pickle(desc_file)
        self.output_name = output_name
        self.df_reactions = self.filter_df_reactions(reactions_file)

        if partitioning_scheme == "hirshfeld":
            self.df_desc = self.df_desc[
                [
                    "smiles",
                    "hirshfeld_charges",
                    "hirshfeld_spin_density",
                    "hirshfeld_charges_fukui_elec",
                    "hirshfeld_charges_fukui_neu",
                    "NMR",
                    "hirshfeld_spin_density_multiplicity",
                    "bond_index_matrix",
                    "Coords",
                ]
            ]

        self.df_desc = self.df_desc.rename(
            columns={
                "hirshfeld_charges": "partial_charge",
                "hirshfeld_spin_density": "spin_dens",
                "hirshfeld_charges_fukui_elec": "fukui_elec",
                "hirshfeld_charges_fukui_neu": "fukui_neu",
                "hirshfeld_spin_density_multiplicity": "spin_dens_triplet",
            }
        )

        if morfeus_file:
            morfeus_df = pd.read_csv(morfeus_file).set_index("smiles")
            morfeus_df["sasa"] = morfeus_df["sasa"].apply(
                lambda x: np.array(list(map(float, str(x).strip("[]\n").split())))
            )
            morfeus_df["pint"] = morfeus_df["pint"].apply(
                lambda x: np.array(list(map(float, str(x).strip("[]\n").split())))
            )
            morfeus_dict = morfeus_df.to_dict()
            self.df_desc["sasa"] = self.df_desc["smiles"].apply(
                lambda x: include_morfeus_data(x, morfeus_dict["sasa"])
            )
            self.df_desc["pint"] = self.df_desc["smiles"].apply(
                lambda x: include_morfeus_data(x, morfeus_dict["pint"])
            )

    def return_valid_df_atoms(self):
        return self.df_reactions

    def filter_df_reactions(self, reactions_file):
        compound_desc_list = self.df_desc["smiles"].values.tolist()
        df_reactions = pd.read_csv(reactions_file)
        df_reactions["desc_avail"] = df_reactions["rxn_smiles"].apply(
            lambda x: verify_smiles(x, compound_desc_list)
        )

        return df_reactions[df_reactions["desc_avail"]]

    def reorder_df_desc(self):
        numbered_smiles_list = []
        for rxn_smiles in self.df_reactions["rxn_smiles"].values.tolist():
            for smiles in rxn_smiles.split(">")[0].split("."):
                numbered_smiles_list.append(smiles)

        desc_dict = self.df_desc.set_index("smiles").T.to_dict()

        df_desc_reordered = pd.DataFrame(
            list(set(numbered_smiles_list)), columns=["smiles"]
        )
        df_desc_reordered["unnumbered_smiles"] = df_desc_reordered["smiles"].apply(
            lambda x: strip_atom_map_num(x)
        )
        df_desc_reordered["ordering"] = df_desc_reordered.apply(
            lambda x: Chem.MolFromSmiles(x["smiles"]).GetSubstructMatch(
                Chem.MolFromSmiles(x["unnumbered_smiles"])
            ),
            axis=1,
        )

        final_columns = ["smiles"]

        for column in self.df_desc.columns:
            if column not in ["smiles", "bond_index_matrix", "Coords"]:
                df_desc_reordered[column] = df_desc_reordered.apply(
                    lambda x: reorder_entry(
                        desc_dict[x["unnumbered_smiles"]][column], x["ordering"]
                    ),
                    axis=1,
                )
                final_columns.append(column)
            elif column == "bond_index_matrix":
                df_desc_reordered["bond_order"] = df_desc_reordered.apply(
                    lambda x: get_bond_orders(
                        x["smiles"],
                        desc_dict[x["unnumbered_smiles"]][column],
                        x["ordering"],
                    ),
                    axis=1,
                )
                final_columns.append("bond_order")
            elif column == "Coords":
                df_desc_reordered["bond_length"] = df_desc_reordered.apply(
                    lambda x: get_bond_lengths(
                        x["smiles"],
                        desc_dict[x["unnumbered_smiles"]][column],
                        x["ordering"],
                    ),
                    axis=1,
                )
                final_columns.append("bond_length")

        self.df_desc = df_desc_reordered[final_columns]
        self.df_desc.dropna(inplace=True)

    def output_atom_descs_chemprop(
        self,
        descriptor_list=[
            "partial_charge",
            "spin_dens",
            "fukui_elec",
            "NMR",
            "fukui_neu",
            "spin_dens_triplet",
        ],
    ):
        df_chemprop = self.df_desc[[name for name in descriptor_list]]
        for name in descriptor_list:
            df_chemprop[name] = turn_elements_into_arrays(df_chemprop[name])

        df_chemprop.to_pickle(f"atom_desc_{self.output_name}_chemprop.pkl", protocol=4)
        df_chemprop.to_csv(f"atom_desc_{self.output_name}_chemprop.csv")

    def output_atom_descs_wln(
        self,
        descriptor_list=[
            "partial_charge",
            "spin_dens",
            "fukui_elec",
            "NMR",
            "fukui_neu",
            "spin_dens_triplet",
            "bond_order",
            "bond_length",
        ],
    ):
        df_wln = self.df_desc[["smiles"] + [name for name in descriptor_list]]

        df_wln.to_pickle(f"atom_desc_{self.output_name}_wln.pkl", protocol=4)
        df_wln.to_csv(f"atom_desc_{self.output_name}_wln.csv")
