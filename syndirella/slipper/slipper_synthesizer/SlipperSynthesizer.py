#!venv/bin/env python3
"""
slipper_synthesizer/CobblersWorkshop.py

This module contains the SlipperSynthesizer class.
"""
from typing import (List, Dict, Tuple)
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem.Fingerprints import FingerprintMols
from syndirella.cobblers_workshop.Library import Library
from syndirella.slipper.slipper_synthesizer.Labeler import Labeler
import pandas as pd
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import os
import shortuuid


class SlipperSynthesizer:
    """
    This class is used to perform the whole process of finding products of the analogues of reactants.
    Since the final elaborated products are 'slippers' in this analogy, the SlipperSynthesizer
    is where these slippers are made.

    This is supposed to be instantiated for each step in the route.
    """

    def __init__(self,
                 library: Library,
                 output_dir: str,
                 atom_ids_expansion: dict = None,
                 additional_info: dict = None):
        self.library = library
        self.output_dir = output_dir
        self.analogues_dataframes_to_react: Dict[str, pd.DataFrame] = {}
        self.analogue_columns: List[str] = None
        self.products: pd.DataFrame = None
        self.reactant_combinations: pd.DataFrame = None
        self.final_products_csv_path: str = None
        self.atom_ids_expansion: dict = atom_ids_expansion
        self.uuid = shortuuid.ShortUUID().random(length=6)
        self.additional_info = additional_info
        self.current_step: int = library.current_step
        self.num_steps: int = library.num_steps

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
        if len(self.analogues_dataframes_to_react) == 1:
            self.products = self.get_products_from_single_reactant()
            return self.products
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
            csv_path = os.path.join(self.output_dir, 'extra', csv_name)
            if os.path.exists(csv_path):
                print(f"Products already exist at {csv_path}. "
                      f"Loading from file...")
                return True
        else:
            final_csv_path = os.path.join(self.output_dir, csv_name)
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
            self.products = pd.read_csv(f"{self.output_dir}/extra/{csv_name}")
        else:
            self.products = pd.read_csv(f"{self.output_dir}/{csv_name}")
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
            #df = self._filter_out_repeats(df, analogue_columns, reactant_prefix)
            df = self.order_analogues(df, reactant_prefix)
            self.analogues_dataframes_to_react[key] = df
        # Filters analogue df by size, shortens if necessary
        self.filter_analogues_by_size()

    def order_analogues(self, df: pd.DataFrame, reactant_prefix: str) -> pd.DataFrame:
        """
        This function is used to order the analogues dataframes by num atom diff to base compound,
        number of reactant matches found, and lead time.
        """
        print(f"Ordering analogues of {reactant_prefix} before finding products...")
        # Add num_atom_diff to base reactant, which is the first reactant
        base_reactant = df[f"{reactant_prefix}_mol"].iloc[0]
        df[f'{reactant_prefix}_num_atom_diff'] = (
            df[f"{reactant_prefix}_mol"].apply(lambda x: self.calc_num_atom_diff_absolute(base_reactant, x)))
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
        df = df[df[orig_r_column]]
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
        if len(self.analogues_dataframes_to_react) < 2:
            if len(list(self.analogues_dataframes_to_react.values())[0]) > 10000:
                print(f"Too many analogues for {list(self.analogues_dataframes_to_react.keys())[0]}.")
                self.analogues_dataframes_to_react[list(self.analogues_dataframes_to_react.keys())[0]] = (
                    self.analogues_dataframes_to_react[list(self.analogues_dataframes_to_react.keys())[0]].head(10000))
            return
        max_allowed_size = 10000
        lengths: List[int] = [len(df) for df in self.analogues_dataframes_to_react.values()]
        product_of_lengths = lengths[0] * lengths[1]
        if product_of_lengths <= max_allowed_size:
            return  # No need to filter
        max_length_each = int(max_allowed_size ** 0.5)  # Taking the square root will give an approximation
        if lengths[0] > max_length_each and lengths[1] <= max_length_each:
            # Cut the first dataframe
            analogue_prefix = list(self.analogues_dataframes_to_react.keys())[0]
            print(f"Too many analogues for {analogue_prefix}.")
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
        if len(self.analogues_dataframes_to_react) < 2:
            return list(self.analogues_dataframes_to_react.values())[0]
        # Get all the analogues dataframes
        r1 = self.analogues_dataframes_to_react['r1']
        r2 = self.analogues_dataframes_to_react['r2']
        combinations = pd.MultiIndex.from_product([r1.index, r2.index], names=['r1', 'r2']).to_frame(index=False)
        # before merging drop analogue_columns
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
        try:
            products = products.explode('combined').reset_index(drop=True)
            # Attempt to split the 'combined' column
            new_columns = pd.DataFrame(products['combined'].tolist(), columns=['smiles', 'num_atom_diff'], index=products.index)
            products[['smiles', 'num_atom_diff']] = new_columns
            products.drop('combined', axis=1, inplace=True)
        except ValueError as e:
            print(f"Error: {e}")
        # Filter products
        products = self.filter_products(products)
        # Add metadata
        products = self.add_metadata(products)
        # Enumerate stereoisomers.
        all_products = self.enumerate_stereoisomers(products)
        print(f"Found {len(set(list(all_products['name'])))} unique products.")
        return all_products

    def get_products_from_single_reactant(self) -> pd.DataFrame:
        """
        This function gets the products from a single reactant (like deprotections).
        """
        reactant = list(self.analogues_dataframes_to_react.values())[0]
        products: pd.DataFrame = reactant.apply(self.apply_reaction_single, axis=1)
        try:
            products = products.explode('combined').reset_index(drop=True)
            # Attempt to split the 'combined' column
            new_columns = pd.DataFrame(products['combined'].tolist(), columns=['smiles', 'num_atom_diff'], index=products.index)
            products[['smiles', 'num_atom_diff']] = new_columns
            products.drop('combined', axis=1, inplace=True)
        except ValueError as e:
            print(f"Error: {e}")
        # Filter products
        products = self.filter_products(products)
        # Add metadata
        products = self.add_metadata(products)
        # Enumerate stereoisomers.
        all_products = self.enumerate_stereoisomers(products)
        print(f"Found {len(set(list(all_products['name'])))} unique products.")
        return all_products

    def apply_reaction_single(self, row) -> pd.Series:
        """
        For mono-molecular reactions:
        This function applies the original reaction to each row of the reactant combinations dataframe. Can return
        multiple products.
        """
        # only get num_atom_diff if final step of route
        reaction: Chem.rdChemReactions = self.library.reaction.reaction_pattern
        r1: str = row['r1_mol']
        products = reaction.RunReactants((r1,))
        if len(products) == 0:
            print("No products found.")
            row['flag'] = None
            row['combined'] = [(None, None)]
            return row
        elif len(products) > 1 or len(products[0]) > 1:
            # check if all products can be sanitized, only keep the ones that can
            row_smiles = []
            row_num_atom_diff = []
            for product in products:
                if self.can_be_sanitized(product[0]):
                    row_smiles.append(Chem.MolToSmiles(product[0]))
                    row_num_atom_diff.append(self.calc_num_atom_diff_absolute(self.library.reaction.product, product[0]))
                if len(row_smiles) > 1:  # if more than 1 product can be sanitized then flag
                    #print(f"Found multiple products at {row.name}. Flagging...")
                    row['flag'] = 'one_of_multiple_products'
                row['combined'] = list(zip(row_smiles, row_num_atom_diff))
                return row
            product = products[0][0]
            if self.can_be_sanitized(product):
                base = self.library.reaction.product
                num_atom_diff = self.calc_num_atom_diff_absolute(base, product)
                row['flag'] = None
                row['combined'] = [(Chem.MolToSmiles(product), num_atom_diff)]
                return row

    def apply_reaction(self, row) -> pd.Series:
        """
        For bimolecular reactions:
        This function applies the original reaction to each row of the reactant combinations dataframe. Can return 
        multiple products. 
        """
        # only get num_atom_diff if final step of route
        reaction: Chem.rdChemReactions = self.library.reaction.reaction_pattern
        r1: str = row['r1_mol']
        r2: str = row['r2_mol']
        products = reaction.RunReactants((r1, r2))
        if len(products) == 0:
            print("No products found.")
            row['flag'] = None
            row['combined'] = [(None, None)]
            # row['smiles'] = None
            # row['num_atom_diff'] = None
            return row
        elif len(products) > 1 or len(
                products[0]) > 1:  # should only return 1 product, if more than 1 then there are selectivity issues
            # check if all products can be sanitized, only keep the ones that can
            row_smiles = []
            row_num_atom_diff = []
            for product in products:
                if self.can_be_sanitized(product[0]):  # only keep products that can be sanitized
                    row_smiles.append(Chem.MolToSmiles(product[0]))
                    row_num_atom_diff.append(self.calc_num_atom_diff_absolute(self.library.reaction.product, product[0]))
            if len(row_smiles) > 1:  # if more than 1 product can be sanitized then flag
                #print(f"Found multiple products at {row.name}. Flagging...")
                row['flag'] = 'one_of_multiple_products'
            row['combined'] = list(zip(row_smiles, row_num_atom_diff))
            return row
        product = products[0][0]
        if self.can_be_sanitized(product):
            base = self.library.reaction.product
            num_atom_diff = self.calc_num_atom_diff_absolute(base, product)
            row['flag'] = None
            row['combined'] = [(Chem.MolToSmiles(product), num_atom_diff)]
            return row

    def can_be_sanitized(self, mol: Chem.Mol):
        assert type(
            mol) == Chem.Mol, f"Expected a Chem.Mol object, got {type(mol)}."  # Make sure it's a Chem.Mol object
        try:
            Chem.SanitizeMol(mol)
            return True
        except:
            return False

    def calc_num_atom_diff_mcs(self, base: Chem.Mol, product: Chem.Mol) -> int:
        """
        This function is used to calculate the number of atoms added to product
        by finding the maximum common substructure (MCS) and then finding the difference in length.
        """
        mcs = rdFMCS.FindMCS([base, product])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        mcs_atoms = mcs_mol.GetNumAtoms()
        new_mol_atoms = product.GetNumAtoms()
        difference = new_mol_atoms - mcs_atoms
        return difference

    def calc_num_atom_diff_absolute(self, base: Chem.Mol, product: Chem.Mol) -> int:
        """
        This function calculates the absolute number of atoms difference between the base and product.
        """
        difference = product.GetNumAtoms() - base.GetNumAtoms()
        return difference

    def filter_products(self, products: pd.DataFrame) -> pd.DataFrame:
        """
        This function is used to filter the products dataframe to remove any rows with None values. Also
        removes duplicates.
        """
        products.dropna(subset=['smiles'], inplace=True, axis=0, how='any')
        products.drop_duplicates(inplace=True, ignore_index=True)
        # reorder by num_atom_diff if calculated
        if 'num_atom_diff' in products.columns:
            products.sort_values(by=['num_atom_diff'], inplace=True)
        products.reset_index(drop=True, inplace=True)
        return products

    def calculate_fingerprints(self, products):
        """Calculate morgan fingerprints for each molecule."""
        products['mol'] = products['smiles'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notnull(x) else None)
        products['fp'] = products['mol'].apply(
            lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024) if x is not None else None)
        return products

    def find_similarity_groups(self, products: pd.DataFrame) -> (pd.DataFrame, int):
        if 'group_id' not in products.columns:
            products['group_id'] = -1  # Initialize all to -1 to indicate no group yet
        fps = products['fp'].tolist()
        n = len(fps)
        # Calculate the fingerprint for the library's reaction product
        product_fp = AllChem.GetMorganFingerprintAsBitVect(self.library.reaction.product, 2, nBits=1024)
        groups = {}
        existing_group_ids = set(products['group_id'])
        group_id = max(existing_group_ids) + 1 if existing_group_ids else 0
        base_group_id = -1  # Initialize with a value to indicate not found
        for i in range(n):
            if fps[i] is None:
                continue
            # Check similarity with the reaction product
            similarity_to_product = DataStructs.FingerprintSimilarity(fps[i], product_fp)
            if similarity_to_product == 1:
                # This is a match; assign it to a group with the reaction product
                if base_group_id == -1:  # If not already assigned
                    while group_id in existing_group_ids:
                        group_id += 1
                    base_group_id = group_id
                    existing_group_ids.add(group_id)
                    group_id += 1
                groups[i] = base_group_id
            for j in range(i + 1, n):
                if fps[j] is None:
                    continue
                similarity = DataStructs.FingerprintSimilarity(fps[i], fps[j])
                if similarity == 1:
                    if i not in groups and j not in groups:
                        while group_id in existing_group_ids:
                            group_id += 1
                        groups[i] = groups[j] = group_id
                        existing_group_ids.add(group_id)
                        group_id += 1
                    elif i in groups and j not in groups:
                        groups[j] = groups[i]
                    elif i not in groups and j in groups:
                        groups[i] = groups[j]
        # Assign group IDs to products based on similarity groups
        for idx, group_id in groups.items():
            products.at[idx, 'group_id'] = group_id
        # For rows without a group, assign a new unique group ID
        ungrouped_indices = products[products['group_id'] == -1].index
        for idx in ungrouped_indices:
            while group_id in existing_group_ids:
                group_id += 1
            products.at[idx, 'group_id'] = group_id
            existing_group_ids.add(group_id)
            group_id += 1
        return products, base_group_id

    def assign_names_based_on_groups(self, products: pd.DataFrame, library_id: str, base_group_id: int) -> pd.DataFrame:
        """Assign names to products based on their group ID, ensuring duplicates have the same name."""
        base_name = f"{library_id}-{self.uuid}-base"
        unique_groups = products['group_id'].unique()
        # Assign names based on group ID
        for group in unique_groups:
            if group == -1:  # Handle molecules without a group
                continue
            if group == base_group_id:
                products.loc[products['group_id'] == group, 'name'] = base_name
                continue
            group_members = products[products['group_id'] == group]
            if not group_members.empty:
                first_name = group_members.iloc[0]['name'] \
                    if pd.notnull(group_members.iloc[0]['name']) \
                    else f"{library_id}-{self.uuid}-{int(group)}"
                products.loc[products['group_id'] == group, 'name'] = first_name
        return products

    def add_metadata(self, products: pd.DataFrame) -> pd.DataFrame:
        print('Adding metadata to products...')
        products['name'] = None  # Add a name column
        products = self.calculate_fingerprints(products)
        products, base_group_id = self.find_similarity_groups(products)
        products = self.assign_names_based_on_groups(products, self.library.id, base_group_id)
        # Add other metadata
        products['reaction'] = self.library.reaction.reaction_name
        products['step'] = self.library.current_step
        products['total_steps'] = self.library.num_steps
        products['base_compound'] = f"{self.library.id}-base"
        products['uuid'] = self.uuid
        # Add additional info
        if self.additional_info:
            for key, value in self.additional_info.items():
                products[key] = value
        products.drop(['mol', 'fp', 'group_id'], axis=1, inplace=True)
        return products

    def enumerate_stereoisomers(self, products: pd.DataFrame) -> pd.DataFrame:
        """
        This function is used to enumerate the stereoisomers of the products.
        """
        # First check if internal step, if yes, don't enumerate stereoisomers
        if self.library.num_steps != self.library.current_step:
            return products
        print("Enumerating stereoisomers since this is the final step...")
        new_rows = []
        for index, row in products.iterrows():
            try:
                stereoisomers = self.find_stereoisomers(row['smiles'])
                for i, iso in enumerate(stereoisomers):
                    new_row = row.copy()
                    new_row['smiles'] = iso
                    new_row['name'] = f"{row['name']}-{chr(65 + i)}"  # Appending A, B, C, etc., to the name
                    new_row['stereoisomer'] = chr(65 + i)
                    new_rows.append(new_row)
            except:
                print(f"Could not enumerate stereoisomers for {row['smiles']}.")
        new_df = pd.DataFrame(new_rows)
        # remove NaNs
        new_df = new_df.dropna(subset=['smiles'])
        new_df.reset_index(drop=True, inplace=True)
        return new_df

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
            print(f"Saving products to {self.output_dir}/extra/{csv_name} \n")
            os.makedirs(f"{self.output_dir}/extra/", exist_ok=True)
            self.products.to_csv(f"{self.output_dir}/extra/{csv_name}", index=False)
        else:
            self.final_products_csv_path: str = f"{self.output_dir}/{csv_name}"
            print(f"Saving final products to {self.final_products_csv_path} \n")
            os.makedirs(f"{self.output_dir}/", exist_ok=True)
            self.products.to_csv(self.final_products_csv_path, index=False)

    def label_products(self):
        """
        This function makes a new instance of the Labeler class and calls the label_products function.
        """
        labeler = Labeler(self.products, self.atom_ids_expansion, self.library)
        self.products = labeler.label_products()
