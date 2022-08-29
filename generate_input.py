from argparse import ArgumentParser, Namespace

from lib import DescFileGenerator


parser = ArgumentParser()
parser.add_argument(
    "--csv_file", type=str, required=True, help=".csv file to read rxn_smiles from"
)
parser.add_argument(
    "--output_name", type=str, required=True, help="base name for output file"
)
parser.add_argument(
    "--prefix", type=str, default="reaction", help="prefix in front of id numbers"
)
parser.add_argument(
    "--side",
    type=str,
    default="reactants",
    help='options: "reactants" (default) or "products"',
)

if __name__ == "__main__":
    args = parser.parse_args()
    file_generator = DescFileGenerator(args.csv_file, args.output_name, args.prefix)
    if args.side == "reactants":
        file_generator.write_reactant_csv()
    if args.side == "products":
        file_generator.write_product_csv()
