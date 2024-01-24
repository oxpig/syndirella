#!venv/bin/env python3
"""
syndirella._cobbler_bench.py

This module contains the CobblerBench class. One instance is made for each step.
"""

import os
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple)
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from ._reaction import Reaction
from ._library import Library
from ..error import ReactionError
from ..smarts import SMARTSHandler

class CobblerBench:
    """
    This class is used to perform the whole process of finding analogues of reactants. Since the final elaborated products
    are 'slippers' in this analogy, the CobblerBench is where these slippers are made.

    It is given a step.
    """
    def __init__(self, product: str, reactants: Tuple[str], reaction_name: str, output_dir: str,
                 smarts_handler: SMARTSHandler, ID: str, num_steps: int, current_step:int):
        self.product: Chem.Mol = Chem.MolFromSmiles(product)
        self.reactants: List[Chem.Mol] = [Chem.MolFromSmiles(reactant) for reactant in reactants]
        self.reaction_name: str = reaction_name
        self.output_dir: str = output_dir
        self.smarts_handler: SMARTSHandler = smarts_handler
        self.ID: str = ID
        self.num_steps: int = num_steps
        self.current_step: int = current_step
        self.reaction: Reaction = None
        self.library = None

    def check_reaction(self):
        """
        This function is used to check the reaction is valid.
        """
        if self.reaction_name not in self.smarts_handler.reaction_smarts.keys():
            raise ReactionError("Reaction name not found in SMARTS handler.")
        if len(self.reactants) != 2:
            raise ReactionError("This function is only for bimolecular reactions.")
        reaction = Reaction(self.product, self.reactants, self.reaction_name, self.smarts_handler)
        return reaction

    def find_analogues_first_step(self):
        """
        This function is used to find the final library of a single step or the first step of a multi step route.
        """
        self.reaction = self.check_reaction()
        # Find attachment ids for all reactants
        self.reaction.find_attachment_ids_for_all_reactants()
        # Find reaction atoms for all reactants
        self.reaction.find_reaction_atoms_for_all_reactants()
        # Find the analogues of reactants
        self.library = Library(self.reaction, self.output_dir, self.ID, self.num_steps, self.current_step)
        self.library.create_library()
        return self.library

    def find_analogues_multi_step(self):
        """
        This function is used to find the analogues of reactants using postera.
        """
        pass


