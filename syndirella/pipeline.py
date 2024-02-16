#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""
import os
from typing import List
import pandas as pd
from rdkit import Chem

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

def _elaborate_compound(smiles: str,
                        compound_set: str,
                        hits: str,
                        template_path: str,
                        hits_path: str,
                        batch_num: int,
                        output_dir: str):
    """
    This function is used to elaborate a single compound.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol, f"Could not create a molecule from the smiles {smiles}."
    # convert hits to a list
    hits = hits.split()
    # create the cobbler
    cobbler = Cobbler(base_compound=smiles, output_dir=output_dir)
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
    print(f"Elaborated compound {smiles}.")
    print()

def run_pipeline(csv_path: str,
                 output_dir: str,
                 template_path: str,
                 hits_path: str,
                 batch_num: int):
    """
    Run the whole syndirella pipeline! 👑
    """
    csv_path = _assert_csv(csv_path)
    df = pd.read_csv(csv_path)
    for index, row in df.iterrows():
        _elaborate_compound(row['smiles'],
                            row['compound_set'],
                            row['hits'],
                            template_path,
                            hits_path,
                            batch_num,
                            output_dir)
    print("Pipeline complete.")



