#!venv/bin/env python3
"""
slipper_synthesizer/_base.py

This module contains the SlipperSynthesizer class.
"""
from typing import (List, Dict, Tuple)
from rdkit import Chem
from rdkit.Chem import rdFMCS
from syndirella.cobblers_workshop._library import Library
from syndirella.slipper.slipper_synthesizer._label import Labeler
import pandas as pd
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import os

class SlipperSynthesizer:
    """
    This class is used to perform the whole process of finding products of the analogues of reactants.
    Since the final elaborated products are 'slippers' in this analogy, the SlipperSynthesizer
    is where these slippers are made.

    This is supposed to be instantiated for each step in the route.
    """

    def __init__(self, library: Library,
                 atom_ids_expansion: dict = None):
        self.library = library
        self.analogues_dataframes_to_react: Dict[str, pd.DataFrame] = {}
        self.analogue_columns: List[str] = None
        self.products: pd.DataFrame = None
        self.reactant_combinations: pd.DataFrame = None
        self.final_products_csv_path: str = None
        self.atom_ids_expansion: dict = atom_ids_expansion

    def get_products(self) -> pd.DataFrame:
        """
        This function is used to find the products of the analogues of reactants. It is the main function that is
        called.
        """
        # Check if products already exist
        if self.check_product_csv_exists():
            self.load_products()
            return self.products
        # Filter analogues
        self.filter_analogues()
        # Get cartesian product of all analogues
        self.reactant_combinations: pd.DataFrame = self.combine_analogues()
        # Find products by applying reaction
        self.products: pd.DataFrame = self.find_products_from_reactants()
        return self.products

    def check_product_csv_exists(self):
        """
        This function checks if the products csv already exists and if so it loads it.
        """
        csv_name = (f"{self.library.id}_{self.library.reaction.reaction_name}_products_"
                    f"{self.library.current_step}of{self.library.num_steps}.csv")
        if self.library.num_steps != self.library.current_step:
            csv_path = os.path.join(self.library.output_dir, 'extra', csv_name)
            if os.path.exists(csv_path):
                print(f"Products already exist at {csv_path}. "
                      f"Loading from file...")
                return True
        else:
            final_csv_path = os.path.join(self.library.output_dir, csv_name)
            if os.path.exists(final_csv_path):
                print(f"Products already exist at {final_csv_path}. "
                      f"Loading from file...")
                return True
        return False

    def load_products(self):
        """
        This function loads the product .csv file.
        """
        csv_name = (f"{self.library.id}_{self.library.reaction.reaction_name}_products_"
                    f"{self.library.current_step}of{self.library.num_steps}.csv")
        if self.library.num_steps != self.library.current_step:
            self.products = pd.read_csv(f"{self.library.output_dir}/extra/{csv_name}")
        else:
            self.products = pd.read_csv(f"{self.library.output_dir}/{csv_name}")
        # Check if any 'Unnamed' columns and remove them
        unnamed_columns = [col for col in self.products.columns if 'Unnamed' in col]
        if len(unnamed_columns) > 0:
            self.products.drop(unnamed_columns, axis=1, inplace=True)
        print(f"Loaded {len(self.products)} products.")

    def filter_analogues(self):
        """
        This function is used to go through the analogue dataframes, passing them to filter_analogues_on_smarts and
        also ordering by metrics.
        Finally it filters the analogues by number, making sure there aren't too many for an obscene number of products.
        """
        for key, value in self.library.analogues_dataframes.items():
            reactant_prefix = key
            df: pd.DataFrame = value[0]
            analogue_columns: Tuple[str, str] = value[1]
            self.analogue_columns = [column for column in analogue_columns]
            df = self.filter_analogues_on_smarts(df, analogue_columns, reactant_prefix)
            df = self.order_analogues(df, reactant_prefix)
            self.analogues_dataframes_to_react[key] = df
        # Filters analogue df by size, shortens if necessary
        self.filter_analogues_by_size()

    def order_analogues(self, df: pd.DataFrame, reactant_prefix: str) -> pd.DataFrame:
        """
        This function is used to order the analogues dataframes by num atom diff to base compound,
        number of reactant matches found, and lead time.
        """
        print(f"Ordering analogues of r{reactant_prefix} before finding products...")
        # Add num_atom_diff to base reactant, which is the first reactant
        base_reactant = df[f"{reactant_prefix}_mol"].iloc[0]
        df[f'{reactant_prefix}_num_atom_diff'] = (
            df[f"{reactant_prefix}_mol"].apply(lambda x: self.calc_num_atom_diff(base_reactant, x)))
        # get columns to sort by
        ordered_columns = [f'{reactant_prefix}_lead_time', f'{reactant_prefix}_num_atom_diff', 'num_matches']
        matching_columns = []
        # Iterate over each substring in the order of preference
        for substring in ordered_columns:
            # Find and append columns that contain the current substring
            matching_columns.extend([col for col in df.columns if substring in col])
        # sort column order
        df = df.sort_values(by=matching_columns, ascending=[True, True, True])
        df.reset_index(drop=True, inplace=True)
        return df

    def filter_analogues_on_smarts(self, df: pd.DataFrame, analogue_columns: Tuple[str, str], reactant_prefix: str) \
            -> pd.DataFrame:
        """
        This function is used to filter the analogues of reactants dataframes to make sure each analogue contains only
        the SMARTS pattern of the original reactant and not the other reactant.
        """
        print('Filtering analogues of reactants on SMARTS...')
        orig_df = df.copy()
        # filter out rows with both 'r1' and 'r2' true (i.e. contains both reactants)
        df = df[~(df[analogue_columns[0]] & df[analogue_columns[1]])]
        # only keep rows with original analogue_prefix true
        orig_r_column = [col for col in analogue_columns if reactant_prefix in col][0]
        df.reset_index(drop=True, inplace=True)
        num_filtered = len(orig_df) - len(df)
        percent_diff = round((num_filtered / len(orig_df)) * 100, 2)
        print(f'Filtered {num_filtered} rows ({percent_diff}%) from {reactant_prefix} dataframe.')
        return df

    def filter_analogues_by_size(self):
        """
        This function is used to filter the analogues dataframes by length. Need to make sure the final combination 
        is less than 10,000.

        If longer than 10,000 will cluster analogues and sample on diversity.
        """
        max_allowed_size = 10000
        lengths: List[int] = [len(df) for df in self.analogues_dataframes_to_react.values()]
        product_of_lengths = lengths[0] * lengths[1]
        if product_of_lengths <= max_allowed_size:
            return  # No need to filter
        max_length_each = int(max_allowed_size ** 0.5)  # Taking the square root will give an approximation
        if lengths[0] > max_length_each and lengths[1] <= max_length_each:
            # Cut the first dataframe
            analogue_prefix = list(self.analogues_dataframes_to_react.keys())[0]
            print(f"Too many analogues for r{analogue_prefix}.")
            analogue_df = self.analogues_dataframes_to_react[analogue_prefix]
            shortened_analogue_df = self.cut_analogues(analogue_df, max_length_each, analogue_prefix)
            self.analogues_dataframes_to_react[analogue_prefix] = shortened_analogue_df
        elif lengths[1] > max_length_each and lengths[0] <= max_length_each:
            # Cut the second dataframe
            analogue_prefix = list(self.analogues_dataframes_to_react.keys())[1]
            print(f"Too many analogues for r{analogue_prefix}.")
            analogue_df = self.analogues_dataframes_to_react[analogue_prefix]
            shortened_analogue_df = self.cut_analogues(analogue_df, max_length_each, analogue_prefix)
            self.analogues_dataframes_to_react[analogue_prefix] = shortened_analogue_df
        else:
            # Cut both dataframes to max_length_each
            print(f"Too many analogues for both reactants.")
            for key in self.analogues_dataframes_to_react.keys():
                self.analogues_dataframes_to_react[key] = self.cut_analogues(
                    self.analogues_dataframes_to_react[key],
                    max_length_each, key)
    def cut_analogues(self, df: pd.DataFrame, max_length_each: int, analogue_prefix: int) -> pd.DataFrame:
        """
        This function is used to cut the analogues dataframes to max_length_each by just taking the head.
        """
        print(f"Cutting {len(df) - max_length_each} analogues from "
              f"{analogue_prefix} dataframe.")
        return df.head(max_length_each)

    def cluster_analogues(self, df: pd.DataFrame, max_length_each: int, analogue_prefix: int) -> pd.DataFrame:
        """
        This function is used to cluster the analogues dataframes to max_length_each by k-means clustering.
        The number of clusters is the number max length each. Might be too much...
        """
        print(f"K-means clustering and sampling {len(df) - max_length_each} analogues from "
              f"r{analogue_prefix} dataframe.")
        # cluster
        df = self.cluster_analogues_on_fingerprint(df)
        # sample
        df = self.sample_analogues(df, max_length_each)
        return df

    def combine_analogues(self):
        """
        This function is used to combine the analogues of reactants into 1 dataframe that the products are found from.
        """
        # Get all the analogues dataframes
        r1 = self.analogues_dataframes_to_react['r1']
        r2 = self.analogues_dataframes_to_react['r2']
        combinations = pd.MultiIndex.from_product([r1.index, r2.index], names=['r1', 'r2']).to_frame(index=False)
        #before merging drop analogue_columns
        r1.drop(self.analogue_columns, axis=1, inplace=True)
        r2.drop(self.analogue_columns, axis=1, inplace=True)
        # merge indicies with original dataframes
        combinations = combinations.merge(r1, left_on='r1', right_index=True)
        combinations = combinations.merge(r2, left_on='r2', right_index=True)
        # drop extra columns
        combinations.drop(['r1', 'r2'], axis=1, inplace=True)
        combinations.reset_index(drop=True, inplace=True)
        # make sure there are no repeats
        combinations.drop_duplicates(inplace=True)
        return combinations

    def find_products_from_reactants(self) -> pd.DataFrame:
        """
        This function is used to find the products of the reactant combinations.
        """
        # Apply reaction to reactant combinations
        products: pd.DataFrame = self.reactant_combinations.apply(self.apply_reaction, axis=1)
        # Explode dataframe in case of multiple products
        products = products.explode('smiles').reset_index(drop=True)
        products = products.explode('num_atom_diff').reset_index(drop=True)
        # Filter products
        products = self.filter_products(products)
        # Add metadata
        products = self.add_metadata(products)
        # Enumerate stereoisomers
        all_products = self.enumerate_stereoisomers(products)
        print(f"Found {len(all_products)} products.")
        return all_products

    def apply_reaction(self, row) -> pd.Series:
        """
        This function applies the original reaction to each row of the reactant combinations dataframe. Can return 
        multiple products. 
        """
        # only get num_atom_diff if final step of route
        reaction: Chem.rdChemReactions = self.library.reaction.reaction_pattern
        r1: str = row['r1_mol']
        r2: str = row['r2_mol']
        products = reaction.RunReactants((r1, r2))
        if len(products) > 1: # should only return 1 product, if more than 1 then there are selectivity issues
            print(f"Found multiple products at {row.name}. Flagging...")
            row['flag'] = 'one_of_multiple_products'
            if all([self.can_be_sanitized(product[0]) for product in products]):
                row['smiles'] = [Chem.MolToSmiles(product[0]) for product in products]
                row['num_atom_diff'] = [self.calc_num_atom_diff(self.library.reaction.product, product[0]) for product in products]
            return row
        if len(products) == 0:
            print("No products found.")
            row['flag'] = None
            row['smiles'] = None
            row['num_atom_diff'] = None
            return row
        product = products[0]
        if self.can_be_sanitized(product):
            base = self.library.reaction.product
            num_atom_diff = self.calc_num_atom_diff(base, product)
            row['flag'] = None
            row['smiles'] = Chem.MolToSmiles(product)
            row['num_atom_diff'] = num_atom_diff
            return row

    def can_be_sanitized(self, mol):
        try:
            Chem.SanitizeMol(mol)
            return True
        except:
            return False

    def calc_num_atom_diff(self, base: Chem.Mol, product: Chem.Mol) -> int:
        """
        This function is used to calculate the number of atoms added to another compound compared to the base compound.
        """
        mcs = rdFMCS.FindMCS([base, product])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        mcs_atoms = mcs_mol.GetNumAtoms()
        new_mol_atoms = product.GetNumAtoms()
        difference = new_mol_atoms - mcs_atoms
        return difference

    def filter_products(self, products: pd.DataFrame) -> pd.DataFrame:
        """
        This function is used to filter the products dataframe to remove any rows with None values. Also
        removes duplicates.
        """
        products.dropna(subset=['smiles'], inplace=True, axis=0)
        products.dropna(inplace=True, axis=1, how='all')
        products.drop_duplicates(inplace=True, ignore_index=True)
        # reorder by num_atom_diff if calculated
        if 'num_atom_diff' in products.columns:
            products.sort_values(by=['num_atom_diff'], inplace=True)
        products.reset_index(drop=True, inplace=True)
        return products

    def add_metadata(self, products: pd.DataFrame) -> pd.DataFrame:
        """
        This function is used to add metadata to the products dataframe.
        """
        products['name'] = products.apply(self.create_name, axis=1)
        products['reaction'] = self.library.reaction.reaction_name
        products['step'] = self.library.current_step
        products['total_steps'] = self.library.num_steps
        products['base_compound'] = f"{self.library.id}_base"
        return products

    def create_name(self, row):
        """
        This function is used to create a name for the product. Does not create a unique name for non final step
        products.
        """
        if row['smiles'] is None or self.library.num_steps != self.library.current_step:
            return None
        else:
            if row['num_atom_diff'] == 0:
                name = f"{self.library.id}-base"
            else:
                name = f"{self.library.id}-{row.name}"
        return name

    def enumerate_stereoisomers(self, products: pd.DataFrame) -> pd.DataFrame:
        """
        This function is used to enumerate the stereoisomers of the products.
        """
        new_rows = []
        for index, row in products.iterrows():
            stereoisomers = self.find_stereoisomers(row['smiles'])
            for i, iso in enumerate(stereoisomers):
                new_row = row.copy()
                new_row['smiles'] = iso
                new_row['name'] = f"{row['name']}-{chr(65 + i)}"  # Appending A, B, C, etc., to the name
                new_row['conformer'] = chr(65 + i)
                new_rows.append(new_row)
        new_df = pd.DataFrame(new_rows)
        exploded_df = pd.concat([products, new_df], ignore_index=True)
        # remove NaNs
        exploded_df = exploded_df.dropna()
        exploded_df.reset_index(drop=True, inplace=True)
        return exploded_df

    def find_stereoisomers(self, smiles: str) -> list():
        # This function should return a list of stereoisomers for the given SMILES string.
        mol = Chem.MolFromSmiles(smiles)
        # Generate all stereoisomers
        opts = StereoEnumerationOptions(unique=True)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        isomer_list = [Chem.MolToSmiles(isomer, isomericSmiles=True) for isomer in isomers]
        return isomer_list

    def save_products(self):
        """
        This function is used to save the products dataframe as a .csv file.
        """
        csv_name = (f"{self.library.id}_{self.library.reaction.reaction_name}_products_"
                    f"{self.library.current_step}of{self.library.num_steps}.csv")
        if self.library.num_steps != self.library.current_step:
            print("Since these products are not the final products they will be saved in the /extra folder. \n")
            print("Saving products to {self.library.output_dir}/extra/{csv_name} \n")
            os.makedirs(f"{self.library.output_dir}/extra/", exist_ok=True)
            self.products.to_csv(f"{self.library.output_dir}/extra/{csv_name}")
        else:
            self.final_products_csv_path: str = f"{self.library.output_dir}/{csv_name}"
            print(f"Saving final products to {self.final_products_csv_path} \n")
            os.makedirs(f"{self.library.output_dir}/", exist_ok=True)
            self.products.to_csv(self.final_products_csv_path)

    def label_products(self):
        """
        This function makes a new instance of the Labeler class and calls the label_products function.
        """
        labeler = Labeler(self.products, self.atom_ids_expansion, self.library)
        self.products = labeler.label_products()











