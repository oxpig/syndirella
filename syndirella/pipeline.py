#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""
import os
from typing import List
import pandas as pd
from rdkit import Chem
import datetime

from syndirella.Cobbler import Cobbler
from syndirella.cobblers_workshop.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper

def _assert_csv(csv_path: str) -> str:
    """
    Make sure that the csv path exists.
    """
    assert os.path.exists(csv_path), "The csv path does not exist."
    # make sure that the csv contains the correct headers
    df = pd.read_csv(csv_path)
    required_columns = ['smiles', 'compound_set', 'hits']
    for col in required_columns:
        assert col in df.columns, f"The csv must contain the column {col}."
    # TODO: could assert that the hits column is a space seperated list of strings

    return csv_path

def _elaborate_compound_with_manual_routes(smiles: str,
                                           hits: str,
                                           template_path: str,):
    """
    This function is used to elaborate a single compound.
    """
    pass


def _elaborate_compound_full_auto(smiles: str,
                        hits: str,
                        template_path: str,
                        hits_path: str,
                        batch_num: int,
                        output_dir: str,
                        additional_info: List[str] = []):
    """
    This function is used to elaborate a single compound.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol, f"Could not create a molecule from the smiles {smiles}."
    # convert hits to a list
    hits = hits.split()
    # get additional info in the correct format
    # create the cobbler
    cobbler = Cobbler(base_compound=smiles,
                      output_dir=output_dir,
                      additional_info=additional_info)
    # get the cobbler workshops
    cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
    for workshop in cobbler_workshops:
        final_library = workshop.get_final_library()
        slipper = Slipper(final_library,
                          template_path,
                          hits_path,
                          hits,
                          batch_num)
        slipper.get_products()
        slipper.place_products()
    print(f"Finished elaborating compound {smiles} at {datetime.datetime.now()}.")
    print()

def run_pipeline(csv_path: str,
                 output_dir: str,
                 template_path: str,
                 hits_path: str,
                 batch_num: int,
                 additional_columns: List[str] = [],
                 manual_routes: bool = False,
                 ):
    """
    Run the whole syndirella pipeline! 👑
    """
    csv_path = _assert_csv(csv_path)
    df = pd.read_csv(csv_path)
    # assert that the df contains additional_columns
    for col in additional_columns:
        assert col in df.columns, f"The csv must contain the column {col}."
    if not manual_routes:
        print("Running the full auto pipeline.")
        for index, row in df.iterrows(): # could make this a parallel for loop
            # TODO: finish this function
            _elaborate_compound_full_auto(row['smiles'],
                                          row['hits'],
                                          template_path,
                                          hits_path,
                                          batch_num,
                                          output_dir,
                                          additional_info=additional_columns)
    else:
        print("Running the pipeline with manual routes.")
        for index, row in df.iterrows():
            _elaborate_compound_with_manual_routes(row['smiles'],
                                                   row['hits'],
                                                   template_path,
                                                   additional_columns)

    print("Pipeline complete.")



