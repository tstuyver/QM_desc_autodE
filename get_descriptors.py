from argparse import ArgumentParser
import pandas as pd

from lib import RxnDescExtractor, AtomDescExtractor


parser = ArgumentParser()
parser.add_argument('--desc_file', type=str, required=True,
                    help='.pkl file to read descs from')
parser.add_argument('--reactions_file', type=str, required=True,
                    help='.csv file to read rxn_smiles from')
parser.add_argument('--output_name', type=str, required=True,
                    help='base name for output file')
parser.add_argument('--morfeus_file', type=str, required=False,
                    help='.csv file to read morfeus descs from')
parser.add_argument('--format', type=str, default='wln',
                    help='options: "wln" (default) or "chemprop"')


if __name__ == '__main__':
    args = parser.parse_args()

    atom_desc_extractor = AtomDescExtractor(args.desc_file, args.reactions_file, args.output_name, args.morfeus_file)
    atom_desc_extractor.reorder_df_desc()

    if args.format == 'wln':
        if args.morfeus_file:
            atom_desc_extractor.output_atom_descs_wln(['partial_charge', 'spin_dens', 'fukui_elec', 'NMR',
                                                            'fukui_neu', 'spin_dens_triplet', 'sasa', 'pint'])
        else:
            atom_desc_extractor.output_atom_descs_wln()

    if args.format == 'chemprop':
        if args.morfeus_file:
            atom_desc_extractor.output_atom_descs_chemprop(['partial_charge', 'spin_dens', 'fukui_elec', 'NMR',
                                                            'fukui_neu', 'spin_dens_triplet', 'sasa', 'pint'])
        else:
            atom_desc_extractor.output_atom_descs_chemprop()

    reaction_desc_extractor = RxnDescExtractor(args.desc_file, args.reactions_file, args.output_name)
    reaction_desc_extractor.determine_reaction_descs()
    
    if args.format == 'wln':
        reaction_desc_extractor.output_reaction_descs_wln()
    if args.format == 'chemprop':
        reaction_desc_extractor.output_reaction_descs_chemprop()
    
    df_reactions_atom = atom_desc_extractor.return_valid_df_atoms()
    df_reactions_react = reaction_desc_extractor.return_valid_df_reactions()
    df_reactions = pd.merge(df_reactions_atom[["rxn_smiles"]], df_reactions_react, how="inner", on=["rxn_smiles"])
    df_reactions.to_csv(f"{args.output_name}_data.csv")
