#!/usr/bin/env python3
"""
syndirella.slipper._base.py

This module contains the Slipper class. A slipper in this metaphor is the set of molecules that is the
product of a reaction.
"""
from typing import (List, Dict, Tuple, Union, Optional)
from rdkit import Chem
import pandas as pd
from syndirella.slipper.slipper_synthesizer._base import SlipperSynthesizer
from syndirella.cobblers_workshop._library import Library

class Slipper:
    """
    This class is instantiated to represent all products for a route.
    """
    def __init__(self, final_library: Library):
        self.products: pd.DataFrame = None
        self.products_smiles: List[str] = None
        self.products_mols: List[Chem.Mol] = None
        self.final_library: Library = final_library

    def get_products(self):
        """
        Main entry to the Slipper class. This function is used to get the products from the final library.
        """
        slipper_synth = SlipperSynthesizer(self.final_library)
        self.products: pd.DataFrame = slipper_synth.get_products()
        slipper_synth.save_products()
        return self.products

    def place_products(self):
        """
        This function is used to place the products with Fragmenstein.
        """
        # add Fragmenstein information to the products


        pass





