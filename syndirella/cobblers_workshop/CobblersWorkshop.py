#!/usr/bin/env python3
"""
syndirella.cobblers_workshop.CobblersWorkshop.py

This module contains the CobblersWorkshop class. One instance of this object is used to describe a full route.
"""

from .CobblerBench import CobblerBench
from typing import (List, Tuple)
from syndirella.cobblers_workshop.Library import Library
from syndirella.SMARTSHandler import SMARTSHandler
from rdkit import Chem
from rdkit.Chem import rdinchi
from syndirella.slipper.Slipper import Slipper
import traceback


class CobblersWorkshop():
    """
    This is the CobblersWorkshop class. It represents a full route.
    """

    def __init__(self,
                 product: str,
                 reactants: List[Tuple],
                 reaction_names: List[str],
                 num_steps: int,
                 output_dir: str,
                 filter: bool,
                 atoms_ids_expansion: dict = None):
        self.product: str = product
        self.id: str = CobblersWorkshop.generate_inchi_ID(self.product)
        self.reactants: List[Tuple[str]] = reactants
        self.reaction_names: List[str] = reaction_names
        self.num_steps: int = num_steps
        self.output_dir: str = output_dir
        self.smarts_handler = SMARTSHandler()
        self.filter: bool = filter
        self.atoms_ids_expansion: dict = atoms_ids_expansion  # should only be internal step
        self.cobbler_benches: List[CobblerBench] = []  # To store instances for each step
        self.first_library: Library = None
        self.final_library: Library = None

    def get_final_library(self):
        """
        This function is used to get the final library of products. It dynamically handles any number of steps.
        """
        try:
            current_library = None
            for step in range(self.num_steps):
                print(f"Step {step + 1} in this route using {self.reaction_names[step]}")
                reactants = self.reactants[step]
                reaction_name = self.reaction_names[step]
                cobbler_bench = CobblerBench(self.product,
                                             reactants,
                                             reaction_name,
                                             self.output_dir,
                                             self.smarts_handler,
                                             self.id,
                                             self.num_steps,
                                             step + 1,
                                             self.filter)
                self.cobbler_benches.append(cobbler_bench)
                current_library = cobbler_bench.find_analogues()
                if step+1 < self.num_steps: # you're not at the last step
                    slipper = Slipper(library=current_library,
                                      atoms_ids_expansion=self.atoms_ids_expansion)
                    slipper.get_products()
                # Update the final library at each step
                self.final_library = current_library
            return self.final_library
        except Exception as e:
            tb = traceback.format_exc()
            print(f"An error occurred in the route elaboration: {e}")
            print(tb)
            return None

    @staticmethod
    def generate_inchi_ID(smiles: str) -> str:
        """
        This function is used to generate a unique id for the route just using the product.
        """
        assert Chem.MolFromSmiles(smiles), f"Could not create a molecule from the smiles {smiles}."
        ID = rdinchi.MolToInchi(Chem.MolFromSmiles(smiles))
        id = rdinchi.InchiToInchiKey(ID[0])
        return id


