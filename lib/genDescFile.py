import pandas as pd

class DescFileGenerator:
    def __init__(self, input_file, output_file, prefix):
        self.df = pd.read_csv(input_file)
        self.output_file = output_file
        self.prefix = prefix

        self.reactant_set = set()
        self.product_set = set()

        self.get_reactant_and_product_set()

    def get_reactant_and_product_set(self):
        self.df["rxn_smiles"].apply(lambda x: self.append_molecule_smiles(x))

    def append_molecule_smiles(self, rxn_smiles):
        reactants, products = rxn_smiles.split(">>")

        for reactant in reactants.split("."):
            self.reactant_set.add(reactant)
        for product in products.split("."):
            self.product_set.add(product)

    def get_dataframe_from_set(self, selected_set):
        df = pd.DataFrame(list(selected_set))
        df["id"] = df.index
        df["id"] = df["id"].apply(lambda x: f'{self.prefix}_{x}')
        df.rename(columns={0: "smiles"}, inplace=True)
        df = df[["id", "smiles"]]

        return df

    def write_reactant_csv(self):
        df_react = self.get_dataframe_from_set(self.reactant_set)
        df_react.to_csv(f'{self.output_file}_reactants.csv')

    def write_product_csv(self):
        df_prod = self.get_dataframe_from_set(self.product_set)
        df_prod.to_csv(f'{self.output_file}_products.csv')
