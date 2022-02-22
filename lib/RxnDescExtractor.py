import pandas as pd
import rdkit.Chem as Chem
import numpy as np

hartree = 627.5094740631


def strip_atom_map_num(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms(): atom.SetAtomMapNum(0)
    stripped_smiles = Chem.MolToSmiles(mol)
    
    return stripped_smiles


def turn_elements_into_arrays(df, column_name):
    df[f'{column_name}'] = df[f'{column_name}'].apply(lambda x: np.array([float(x)]))


class RxnDescExtractor:
    def __init__(self, desc_file, reactions_file, output_name):
        self.df_desc = pd.read_pickle(desc_file)
        self.df_reactions = pd.read_csv(reactions_file)
        self.output_name = output_name

        self.mol_desc = self.extract_descriptor_values()

        self.reaction_desc = self.determine_reaction_descs()

    def determine_reaction_descs(self):
        self.df_reactions['G_orb'] = self.df_reactions['rxn_smiles'].apply(lambda x: self.get_promotion_gap_orbitals(x))
        self.df_reactions['G'] = self.df_reactions['rxn_smiles'].apply(lambda x: self.get_promotion_gap_scf(x))
        self.df_reactions['G_alt1_orb'] = self.df_reactions['rxn_smiles'].apply(lambda x: self.get_promotion_gap_alt_orbitals(x, 1))
        self.df_reactions['G_alt1'] = self.df_reactions['rxn_smiles'].apply(lambda x: self.get_promotion_gap_alt_scf(x, 1))
        self.df_reactions['G_alt2_orb'] = self.df_reactions['rxn_smiles'].apply(lambda x: self.get_promotion_gap_alt_orbitals(x, 2))
        self.df_reactions['G_alt2'] = self.df_reactions['rxn_smiles'].apply(lambda x: self.get_promotion_gap_alt_scf(x, 2))
        
    def extract_descriptor_values(self):
        mol_desc = {'reactants_singlet': dict(zip(self.df_desc['smiles'], self.df_desc['SCF'])),
        'reactants_homo': dict(zip(self.df_desc['smiles'], self.df_desc['homo'])),
        'reactants_lumo': dict(zip(self.df_desc['smiles'], self.df_desc['lumo'])),
        'reactants_triplet': dict(zip(self.df_desc['smiles'], self.df_desc['SCF_multiplicity'])),
        'reactants_plus': dict(zip(self.df_desc['smiles'], self.df_desc['SCF_plus1'])),
        'reactants_minus': dict(zip(self.df_desc['smiles'], self.df_desc['SCF_minus1']))}

        return mol_desc

    def get_promotion_gap_orbitals(self, rxn_smiles):
        reactants = rxn_smiles.split(">")[0].split(".")
        promotion_energy = 0
        for reactant in reactants:
            try:
                stripped_reactant = strip_atom_map_num(reactant)
                promotion_energy += (self.mol_desc['reactants_homo'][stripped_reactant] - \
                                    self.mol_desc['reactants_lumo'][stripped_reactant]) * hartree
            except Exception as e:
                print(f'Error for {rxn_smiles}: {e}')
                promotion_energy = None
    
        return promotion_energy

    def get_promotion_gap_scf(self, rxn_smiles):
        reactants = rxn_smiles.split(">")[0].split(".")
        promotion_energy = 0
        for reactant in reactants:
            try:
                stripped_reactant = strip_atom_map_num(reactant)
                promotion_energy += (self.mol_desc['reactants_triplet'][stripped_reactant] - \
                                    self.mol_desc['reactants_singlet'][stripped_reactant]) * hartree
            except Exception as e:
                print(f'Error for {rxn_smiles}: {e}')
                promotion_energy = None
    
        return promotion_energy

    def get_promotion_gap_alt_orbitals(self, rxn_smiles, num):
        reactants = rxn_smiles.split(">")[0].split(".")
        for reactant in reactants:
            reactant_mol = Chem.MolFromSmiles(reactant)
            try:
                if num == 1:
                    if any([atom.GetFormalCharge() == 1 for atom in reactant_mol.GetAtoms()]):
                        dipole = strip_atom_map_num(reactant)
                        electron_affinity = - self.mol_desc['reactants_lumo'][dipole] * hartree
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        ionisation_potential = - self.mol_desc['reactants_homo'][dipolarophile] * hartree
                if num == 2:
                    if any([atom.GetFormalCharge() == 1 for atom in reactant_mol.GetAtoms()]):
                        dipole = strip_atom_map_num(reactant)
                        ionisation_potential = - self.mol_desc['reactants_homo'][dipole] * hartree
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        electron_affinity = - self.mol_desc['reactants_lumo'][dipolarophile] * hartree 
            except Exception as e:
                print(f'error for {rxn_smiles}: {e}')
                return None

        return ionisation_potential - electron_affinity

    def get_promotion_gap_alt_scf(self, rxn_smiles, num):
        reactants = rxn_smiles.split(">")[0].split(".")
        for reactant in reactants:
            reactant_mol = Chem.MolFromSmiles(reactant)
            try:
                if num == 1:
                    if any([atom.GetFormalCharge() == 1 for atom in reactant_mol.GetAtoms()]):
                        dipole = strip_atom_map_num(reactant)
                        electron_affinity = (self.mol_desc['reactants_singlet'][dipole] - 
                                            self.mol_desc['reactants_minus'][dipole]) * hartree
                        if electron_affinity < 0:
                            electron_affinity = \
                                (- (self.mol_desc['reactants_lumo'][dipole] + self.mol_desc['reactants_homo'][dipole]) \
                                - (self.mol_desc['reactants_plus'][dipole] - self.mol_desc['reactants_singlet'][dipole])) * hartree
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        ionisation_potential = (self.mol_desc['reactants_plus'][dipolarophile] - 
                                                self.mol_desc['reactants_singlet'][dipolarophile]) * hartree
                if num == 2:
                    if any([atom.GetFormalCharge() == 1 for atom in reactant_mol.GetAtoms()]):
                        dipole = strip_atom_map_num(reactant)
                        ionisation_potential = (self.mol_desc['reactants_plus'][dipole] - 
                                                self.mol_desc['reactants_singlet'][dipole]) * hartree
                    else:
                        dipolarophile = strip_atom_map_num(reactant)
                        electron_affinity = (self.mol_desc['reactants_singlet'][dipolarophile] - 
                                            self.mol_desc['reactants_minus'][dipolarophile]) * hartree
                        if electron_affinity < 0:
                            electron_affinity = \
                                (- (self.mol_desc['reactants_lumo'][dipolarophile] + self.mol_desc['reactants_homo'][dipolarophile]) \
                                - (self.mol_desc['reactants_plus'][dipolarophile] - self.mol_desc['reactants_singlet'][dipolarophile])) * hartree
            except Exception as e:
                print(f'error for {rxn_smiles}: {e}')
                return None

        return ionisation_potential - electron_affinity 

    def output_reaction_descs_chemprop(self, descriptor_list=['E_r', 'G', 'G_alt1', 'G_alt2']):
        df_chemprop = self.df_reactions[[name for name in descriptor_list]]
        for name in descriptor_list:
            df_chemprop = turn_elements_into_arrays(df_chemprop, name)

        df_chemprop.to_pickle(f'reaction_desc_{self.output_name}_chemprop.pkl')
        df_chemprop.to_csv(f'reaction_desc_{self.output_name}_chemprop.csv')

    def output_reaction_descs_wln(self, descriptor_list=['E_r', 'G', 'G_alt1', 'G_alt2']):
        df_wln = self.df_reactions[['rxn_smiles'] + [name for name in descriptor_list]]
        df_wln = df_wln.rename(columns={'rxn_smiles':'smiles'})

        df_wln.to_pickle(f'reaction_desc_{self.output_name}_wln.pkl')
        df_wln.to_csv(f'reaction_desc_{self.output_name}_wln.csv')
