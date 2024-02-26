#!venv/bin/env python3
"""
syndirella.CobblerBench.py

This module contains the CobblerBench class. One instance is made for each step.
"""

import os
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple)
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from .Reaction import Reaction
from .Library import Library
from ..error import ReactionError
from ..SMARTSHandler import SMARTSHandler

class CobblerBench:
    """
    This class is used to perform the whole process of finding analogues of reactants. Since the final elaborated products
    are 'slippers' in this analogy, the CobblerBench is where these slippers are made.

    It is given a step.
    """
    def __init__(self,
                 product: str,
                 reactants: Tuple[str],
                 reaction_name: str,
                 output_dir: str,
                 smarts_handler: SMARTSHandler,
                 id: str,
                 num_steps: int,
                 current_step: int,
                 filter: bool):
        self.product: Chem.Mol = Chem.MolFromSmiles(product)
        self.reactants: List[Chem.Mol] = [Chem.MolFromSmiles(reactant) for reactant in reactants]
        self.reaction_name: str = reaction_name
        self.output_dir: str = output_dir
        self.smarts_handler: SMARTSHandler = smarts_handler
        self.id: str = id
        self.num_steps: int = num_steps
        self.current_step: int = current_step
        self.filter: bool = filter
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
        reaction = Reaction(self.product,
                            self.reactants,
                            self.reaction_name,
                            self.smarts_handler)
        return reaction

    def find_analogues(self):
        """
        This function is used to find analogues of any step along the route.
        """
        self.reaction = self.check_reaction()
        # Find attachment ids for all reactants
        self.reaction.find_attachment_ids_for_all_reactants()
        # Find reaction atoms for all reactants
        self.reaction.find_reaction_atoms_for_all_reactants()
        # Find the analogues of reactants
        self.library = Library(self.reaction,
                               self.output_dir,
                               self.id,
                               self.num_steps,
                               self.current_step,
                               self.filter)
        self.library.create_library()
        return self.library
