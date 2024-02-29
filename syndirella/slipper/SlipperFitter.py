#!venv/bin/env python3
"""
slipper_fitter/CobblersWorkshop.py

This module contains the SlipperFitter class.
"""
from typing import (List, Dict, Tuple, Union, Optional)
from fragmenstein import Laboratory, Wictor
import pandas as pd
from fragmenstein.laboratory.validator import place_input_validator
from rdkit import Chem
import os, logging
import time


class SlipperFitter:
    """
    This class is instantiated to place all the products in the template using hits_path
    """

    def __init__(self,
                 final_products: pd.DataFrame,
                 template_path: str,
                 hits_path: str,
                 hits_names: List[str],
                 n_cores: int,
                 timeout: int,
                 batch_num: int,
                 output_dir: str,
                 final_products_csv_path: str):
        self.final_products: pd.DataFrame = final_products
        self.template_path: str = template_path
        self.hits_path: str = hits_path
        self.hits_names: List[str] = hits_names
        self.input_df: pd.DataFrame = None
        self.placements: pd.DataFrame = None
        self.merged_placements: pd.DataFrame = None

        # Placements variables set
        self.n_cores: int = n_cores
        self.timeout: int = timeout
        self.batch_num: int = batch_num
        self.output_dir: str = output_dir

        self.final_products_csv_path: str = final_products_csv_path
        self.output_path = None

    def fit_products(self):
        """
        Main entry to the SlipperFitter class. This function is used to fit the products to the final library.
        """
        self.input_df = self.prep_products()
        self.placements = self.place_products()
        self.edit_placements()
        self.merged_placements = self.merge_placements()
        return self.merged_placements

    def prep_products(self) -> pd.DataFrame:
        """
        This function is used to prepare the inputs for Fragmenstein.
        """
        # input_df is a copy of final_products but removing duplicates of names and other metadata from synthesizer step
        input_df: pd.DataFrame = self._prepare_input_df()
        input_df = self.add_hits(input_df)
        # add typing to columns
        input_df = input_df.astype({'name': str, 'smiles': str})
        # Check if there are any duplicates in the input_df
        if input_df.duplicated(subset='name').any():
            raise ValueError('There are duplicates of names in the product dataframe to place.')
        return input_df

    def _prepare_input_df(self) -> pd.DataFrame:
        """
        Creates the input dataframe for the placements. Will remove duplicates of names and other metadata from
        synthesizer step.
        """
        input_df: pd.DataFrame = self.final_products.copy(deep=True)
        # drop duplicates of names
        input_df = input_df.drop_duplicates(subset='name')
        # drop columns that are not needed
        input_df = input_df[['name', 'smiles']]
        # place number in batch
        if self.batch_num > 0:
            print(f'Only placing the top {self.batch_num} products.')
            input_df = input_df.iloc[:self.batch_num]
        self._print_diff(self.final_products, input_df)
        return input_df

    def _print_diff(self, orig_df: pd.DataFrame, input_df: pd.DataFrame):
        """
        This function is used to print the difference between the original number of analogues and the number of
        valid analogues.
        """
        assert len(input_df) <= len(orig_df), ("Problem with finding unique analogues. There are more than were in "
                                               "the original list of analogues")
        num_placed = len(input_df)
        percent = round(((len(input_df) / len(orig_df)) * 100), 2)
        print(f'Placing {len(input_df)} ({percent}%) unique analogues out of {len(orig_df)} analogues.')

    def add_hits(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """
        This function adds the hits_path as mol objects to input_df['hits_path'].
        """
        # load hits_path either from mol or sdf
        if os.path.splitext(self.hits_path)[1] == '.mol':
            print('This is a mol file')
            hits: List[Chem.Mol] = [Chem.MolFromMolFile(self.hits_path.strip())]
        else:
            with Chem.SDMolSupplier(self.hits_path.strip()) as sd:
                hits: List[Chem.Mol] = list(sd)
        # Find which hits are in the hit_names
        hits = [
            hit for hit in hits
            if any(hit_name in hit.GetProp('_Name') for hit_name in self.hits_names)
        ]
        # add hits column with contents of hits_names as space separated string
        input_df['hits'] = input_df.apply(lambda row: hits, axis=1)
        return input_df

    def place_products(self) -> pd.DataFrame:
        """
        This function places products using Fragmenstein.
        """
        start_time = time.time()  # Start timing
        # set up Wictor
        self.setup_Fragmenstein()
        # Get pdbblock from template_path
        with open(self.template_path) as fh:
            pdbblock: str = fh.read()
        lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=False)
        placements: pd.DataFrame = lab.place(place_input_validator(self.input_df), n_cores=self.n_cores,
                                             timeout=self.timeout)
        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time  # Calculate elapsed time
        print(f"Placing {len(self.input_df)} run time: {elapsed_time:.2f} seconds")  # Print the elapsed time
        return placements

    def setup_Fragmenstein(self):
        """
        This function sets up Fragmenstein to run.
        """
        # Using Wictor to place just by RDKit minimisation
        # make output directory
        self.output_path: str = os.path.join(self.output_dir, 'output')
        os.makedirs(self.output_path, exist_ok=True)
        Wictor.work_path = self.output_path
        os.chdir(self.output_dir)  # this does the trick
        Wictor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
        Wictor.monster_joining_cutoff = 5  # Å
        Wictor.quick_reanimation = False  # for the impatient
        Wictor.error_to_catch = Exception  # stop the whole laboratory otherwise
        Wictor.enable_stdout(logging.CRITICAL)
        # Wictor.enable_logfile(os.path.join(self.output_path, f'fragmenstein.log'), logging.ERROR)
        Laboratory.Victor = Wictor

    def edit_placements(self):
        """
        This function edits the placements.
        """
        self.placements.loc[(self.placements['∆∆G'] < 0) & (self.placements['comRMSD'] < 2), 'outcome'] = 'acceptable'
        self.placements.loc[(self.placements['∆∆G'] > -1) &
                            (self.placements.outcome == 'acceptable'), 'outcome'] = 'weak'
        self.placements['unminimized_mol'] = self.placements.unminimized_mol.fillna(Chem.Mol())
        self.fix_intxns()
        if self.placements.outcome.value_counts().get('acceptable') is None:
            num_success = 0
            percent_success = 0
        else:
            num_success = self.placements.outcome.value_counts()['acceptable']
            percent_success = round((self.placements.outcome.value_counts()['acceptable'] / len(self.placements) * 100)
                                    , 2)
        print(f'{num_success} ({percent_success}%) successful placements '
              f'where ∆∆G < -1 and RMSD < 2 Å.')

    def fix_intxns(self):
        intxn_names: List[tuple] = [c for c in self.placements.columns if isinstance(c, tuple)]
        for intxn_name in intxn_names:
            self.placements[intxn_name] = self.placements[intxn_name].fillna(0).astype(int)

    def merge_placements(self) -> pd.DataFrame:
        """
        This function merges the fragmenstein output with products csv.
        """
        # Define the columns to drop
        columns_to_drop = ['smiles', 'binary_hits', 'mode', 'runtime', 'disregarded', 'unmin_binary', 'min_binary',
                           'hit_binaries', 'unminimized_mol', 'minimized_mol', 'hit_mols', 'hit_names']
        # Filter out the columns that do not exist in merged_placements
        columns_to_drop = [col for col in columns_to_drop if col in self.placements.columns]
        # Drop the columns if they exist
        self.placements = self.placements.drop(columns=columns_to_drop)
        # merge on name
        merged_placements = pd.merge(self.final_products, self.placements, on='name', how='left')
        return merged_placements

    def save_placements(self):
        """
        This function saves the placements as a .pkl.gz and .csv file.
        """
        self.placements.to_pickle(os.path.join(self.output_dir, 'fragmenstein_placements.pkl.gz'))
        self.placements.to_csv(os.path.join(self.output_dir, 'fragmenstein_placements.csv'))
        # get final name same as final products csv
        merged_placements_path = self.final_products_csv_path.split('.')[0] + '_placements.csv'
        self.merged_placements.to_csv(merged_placements_path)
