#!venv/bin/env python3
"""
syndirella.cobblers_workshop.Library.py

This module contains the Library class. This class contains information about the analogue library. It will create the
analogue library from the Reaction object. It will also store the analogue library as a .csv file.
"""

import pandas as pd
from rdkit.Chem.FilterCatalog import *

from typing import (List, Dict, Tuple)
import os
from rdkit import DataStructs
from rdkit.Chem import AllChem
import glob
import pickle

from .Reaction import Reaction
from syndirella.Postera import Postera
import syndirella
from syndirella.Fairy import Fairy


class Library:
    """
    This class contains information about the analogue library. It will create the analogue library from the Reaction
    object. It will also store the analogue library as a .csv file.

    """
    def __init__(self, reaction: Reaction, output_dir: str, id: str, num_steps: int, current_step: int, filter: bool):
        self.reaction: Reaction = reaction
        self.id: str = id
        self.output_dir: str = f"{output_dir}/{self.id}"
        self.extra_dir_path: str = f"{self.output_dir}/extra"
        self.num_steps: int = num_steps
        self.current_step: int = current_step
        self.analogues_dataframes: Dict[str: pd.DataFrame] = {}
        self.filter: bool = filter

        self.r1 = None
        self.r2 = None

    def create_library(self):
        """
        This function is used to create the analogue library from the Reaction object.
        """
        for key, value in self.reaction.matched_smarts_to_reactant.items():
            reactant = value[0]
            reactant_smarts = key
            analogue_prefix = value[2]
            df, analogue_columns = self.process_reactant(reactant, reactant_smarts, analogue_prefix)
            df: pd.DataFrame
            analogue_columns: Tuple[str, str]
            self.analogues_dataframes[analogue_prefix] = (df, analogue_columns)

    def process_reactant(self, reactant: Chem.Mol, reactant_smarts: str, analogue_prefix: str) -> pd.DataFrame:
        # search for analogues csv if already created. Perform for all reactants.
        reactant_analogues_path, internal_step = self.check_analogues_csv_exists(analogue_prefix, reactant)
        if reactant_analogues_path is not None:
            analogues_full: Dict[str, float] = self.load_library(reactant_analogues_path, analogue_prefix,
                                                                 internal_step)
            analogues: List[str] = list(analogues_full.keys())
        else: # perform database search
            postera_search = Postera(self.reaction, self.output_dir, search_type="superstructure")
            # returns a dictionary of analogues where values = min lead time found
            analogues_full: Dict[str, float] = postera_search.perform_database_search(reactant)
            analogues: List[str] = list(analogues_full.keys())
        processed_analogues_df, analogue_columns = (
            self.process_analogues(analogues, reactant_smarts, analogue_prefix, analogues_full))
        processed_analogues_df: pd.DataFrame
        analogue_columns: Tuple[str, str]
        self.save_library(processed_analogues_df, analogue_prefix)
        return processed_analogues_df, analogue_columns

    def process_analogues(self, analogues: List[str], reactant_smarts: str, analogue_prefix: str,
                          analogues_full: Dict[str, float] = None) -> pd.DataFrame:
        """
        This function puts list of analogues in dataframe and does SMART checking to check if the analogues contains
        the SMARTS pattern of the original reactant and against all other reactants SMARTS.
        """
        analogues_mols: List[Chem.Mol] = [Chem.MolFromSmiles(analogue) for analogue in analogues]
        analogues_mols: List[Chem.Mol] = Fairy.remove_repeat_mols(analogues_mols)
        self.print_diff(analogues, analogues_mols, analogue_prefix)
        if self.filter:
            analogues: List[str] = self.filter_analogues(analogues, analogue_prefix)
        reactant_smarts_mol: Chem.Mol = Chem.MolFromSmarts(reactant_smarts)

        contains_smarts_pattern, num_matches = self.check_analogue_contains_smarts_pattern(analogues_mols,
                                                                                           reactant_smarts_mol)
        contains_smarts_pattern: List[bool]
        num_matches: List[int]

        contains_other_reactant_smarts_pattern, other_reactant_prefix = (
            self.check_analogue_contains_other_reactant_smarts_pattern(analogues_mols, reactant_smarts))
        contains_other_reactant_smarts_pattern: List[bool]
        other_reactant_prefix: str

        # get lead time from analogues
        analogues = [Chem.MolToSmiles(analogue) for analogue in analogues_mols]
        if analogues_full is not None:
            lead_times = [analogues_full[analogue] for analogue in analogues]
            assert len(lead_times) == len(analogues), "Problem with finding lead times."

        analogues_df = (
            pd.DataFrame({f"{analogue_prefix}_smiles": analogues,
                          f"{analogue_prefix}_mol": analogues_mols,
                          f"{analogue_prefix}_{self.reaction.reaction_name}": contains_smarts_pattern,
                          f"{analogue_prefix}_{self.reaction.reaction_name}_num_matches": num_matches,
                          f"{other_reactant_prefix}_{self.reaction.reaction_name}": contains_other_reactant_smarts_pattern,
                          f"{analogue_prefix}_lead_time": lead_times}))

        if self.filter:
            analogues_df['is_PAINS_A'] = False
        return analogues_df, (f"{analogue_prefix}_{self.reaction.reaction_name}",
                              f"{other_reactant_prefix}_{self.reaction.reaction_name}")

    def check_analogue_contains_smarts_pattern(self, analogues_mols: List[Chem.Mol], reactant_smarts_mol: Chem.Mol):
        """
        This function is used to check if the analogues contains the original reactant. Will return
        dictionary of boolean values and the number of matches.
        """
        print('Checking if analogues contain SMARTS pattern of original reactant...')
        matching = [bool(analogue.GetSubstructMatches(reactant_smarts_mol)) for analogue in analogues_mols]
        num = [len(analogue.GetSubstructMatches(reactant_smarts_mol)) for analogue in analogues_mols]
        assert len(matching) == len(analogues_mols), "Problem with finding matches."
        assert len(num) == len(analogues_mols), "Problem with finding number of matches."
        return matching, num

    def filter_analogues(self, analogues: List[str], analogue_prefix: str) -> List[str]:
        """
        This function is used to filter out analogues.
        """
        mols: List[Chem.Mol] = [Chem.MolFromSmiles(mol) for mol in analogues]
        passing_mols: List[Chem.Mol] = self.filter_on_substructure_filters(mols)
        self.print_diff(mols, passing_mols, analogue_prefix)
        filtered_analogues: List[str] = [Chem.MolToSmiles(mol) for mol in passing_mols]
        return filtered_analogues

    def print_diff(self, mols: List[Chem.Mol], valid_mols: List[Chem.Mol], analogue_prefix: str):
        """
        This function is used to print the difference between the original number of analogues and the number of
        valid analogues.
        """
        assert len(valid_mols) <= len(mols), ("Problem with finding valid molecules. There are more than were in "
                                              "the original list of molecules.")
        num_filtered = len(mols) - len(valid_mols)
        percent_diff = round((num_filtered / len(mols)) * 100, 2)
        print(f'Removed {num_filtered} invalid or repeated molecules ({percent_diff}%) of {analogue_prefix} analogues.')

    def filter_on_substructure_filters(self, mols: List[Chem.Mol], ) -> List[str]:
        """
        This function is used to filter out analogues that do not pass the substructure filters.
        """
        print('Filtering analogues on PAINS_A filters...')
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        catalog = FilterCatalog(params)
        passing_molecules: List[Chem.Mol] = []
        for mol in mols:
            if not catalog.HasMatch(mol):
                passing_molecules.append(mol)
        return passing_molecules

    def check_analogue_contains_other_reactant_smarts_pattern(self, analogues_mols: List[Chem.Mol],
                                                              reactant_smarts: str) -> List[bool]:
        """
        This function is used to check if the analogues contains the SMARTS patterns of the other reactant.
        """
        # get other reactant SMARTS pattern
        print('Checking if analogues contain SMARTS pattern of other reactant...')
        other_reactant_smarts_pattern = [smarts for smarts in self.reaction.matched_smarts_to_reactant.keys() if not smarts == reactant_smarts][0]
        assert other_reactant_smarts_pattern is not None, "Other reactant SMARTS pattern not found."
        other_reactant_prefix = self.reaction.matched_smarts_to_reactant[other_reactant_smarts_pattern][2]
        other_reactant_smarts_mol: Chem.Mol = Chem.MolFromSmarts(other_reactant_smarts_pattern)
        return [bool(analogue.GetSubstructMatches(other_reactant_smarts_mol)) for analogue in analogues_mols], other_reactant_prefix

    def check_analogue_contains_all_smarts_patterns(self, analogues_mols: List[Chem.Mol]) -> List[bool]:
        """
        This function is used to check if the analogues contains all the SMARTS patterns of the other reactants.
        """
        return NotImplementedError()

    def check_analogues_csv_exists(self, analogue_prefix: str, reactant: Chem.Mol) -> str and bool:
        """
        This function is used to check if the analogue library has already been created and saved as a .csv file.
        """
        internal_step = False
        if self.current_step != 1 and self.current_step == self.num_steps:
            print('Since this is the final step, looking for the products .csv from first step...')
            csv_path: str = self.find_products_csv_name(reactant)
            internal_step = True
            return csv_path, internal_step
        os.makedirs(self.extra_dir_path, exist_ok=True)
        csv_name = (f"{self.id}_{self.reaction.reaction_name}_"
                    f"{analogue_prefix}_{self.current_step}of{self.num_steps}.csv")
        if os.path.exists(f"{self.extra_dir_path}/{csv_name}"):
            return f"{self.extra_dir_path}/{csv_name}", internal_step
        else:
            return None, None

    def find_products_csv_name(self, reactant: Chem.Mol) -> str:
        """
        This function is used to find the name of the products .csv file by comparing the reactant to the product with
        bit vector similarity.
        """
        product_csvs: List[str] = glob.glob(f"{self.extra_dir_path}/"
                                            f"{self.id}_*products_{self.current_step-1}of"
                                            f"{self.num_steps}.csv")
        for i, path in enumerate(product_csvs):
            # Find if the reactant is the same as the product
            df = pd.read_csv(path)
            # Iterate through the top 100 products
            for product in df["smiles"][:100]:
                product_mol = Chem.MolFromSmiles(product)
                # Calculate similarity score
                similarity_score = DataStructs.FingerprintSimilarity(
                    AllChem.GetMorganFingerprintAsBitVect(reactant, 2),
                    AllChem.GetMorganFingerprintAsBitVect(product_mol, 2))
                # If a perfect match is found, return the path
                if similarity_score == 1:
                    print(f"Found {path} as the products .csv from first step.")
                    return path
        return None

    def save_library(self, df: pd.DataFrame, analogue_prefix: str):
        """
        This function is used to save the analogue library as a .csv file in self.extra_dir_path
        """
        os.makedirs(f"{self.extra_dir_path}", exist_ok=True)
        csv_name = f"{self.id}_{self.reaction.reaction_name}_{analogue_prefix}_{self.current_step}of{self.num_steps}.csv"
        print(f"Saving {analogue_prefix} analogue library to {self.extra_dir_path}/{csv_name} \n")
        df.to_csv(f"{self.extra_dir_path}/{csv_name}", index=False)

    def load_library(self, reactant_analogues_path: str, analogue_prefix: str, internal_step: bool) -> Dict[str, float]:
        """
        This function is used to load the analogue library from a .csv file. Returns a dictionary of the smiles as keys
        and the lead time (if found) as float. Otherwise return lead time as None.
        """
        try:
            # find column with analogue smiles
            df = pd.read_csv(reactant_analogues_path)
            if not internal_step:
                analogues = df[[f"{analogue_prefix}_smiles", f"{analogue_prefix}_lead_time"]]
                analogues_full = {row[0]: row[1] for row in analogues.itertuples(index=False)}
                return analogues_full
            if internal_step:
                analogues = df["smiles"].tolist()
                analogues_full = {analog: None for analog in analogues}
                return analogues_full
        except KeyError:
            print(f"Could not find analogue column in already existing product.csv at {reactant_analogues_path}. "
                  f"Stopping...")

    def save(self):
        """
        This function saves the library object in the self.extra_output_dir_path.
        """
        with open(f"{self.extra_dir_path}/{self.id}_library.pkl", "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(output_dir: str, product: str) -> 'Library':
        """
        This function loads the library object. Will have to load all library objects
        """
        # TODO: Change this to load all library objects, but check that it contains the specific variable you want
        # get all .pkl in output_dir.
        library_pkls: List[str] = glob.glob(f"{output_dir}/*library*.pkl")
        # get id from product SMILES
        id = syndirella.cobblers_workshop.CobblersWorkshop.generate_inchi_ID(product)
        with open(library_pkls[0], "rb") as f:
            library = pickle.load(f)
            if library.id == id and library.num_steps == library.current_step: # return the final library
                print(f"Loaded {library.id} final library.")
                return library
        # if not found, return None
        print(f"Could not load {id} final library.")
        return None

