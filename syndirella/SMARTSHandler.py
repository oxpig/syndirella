#!/usr/bin/env python3
"""
SMARTSHandler.py

This module contains the SMARTSHandler class. This class contains information about the reaction SMARTS.
"""
import json
from collections import OrderedDict
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple)
from .cli_defaults import cli_default_settings

class SMARTSHandler:
    def __init__(self):
        with open(cli_default_settings['rxn_smarts_path']) as f:
            reaction_smarts = json.load(f)
        self.reaction_smarts = {name: ReactionFromSmarts(val) for name, val in reaction_smarts.items()}
        self.reactant1_dict = OrderedDict()
        self.reactant2_dict = OrderedDict()
        self.product_smarts = OrderedDict()
        for name, smart in reaction_smarts.items():
            reactants, prod = smart.split(">>")
            try:
                react1, react2 = reactants.split(".")
            except ValueError:
                react1, react2 = reactants, None
            self.reactant1_dict[name] = react1
            self.reactant2_dict[name] = react2
            self.product_smarts[name] = prod
        self.pattern_products = self.from_SMARTS_to_patterns(self.product_smarts)
        self.pattern_reactant1 = self.from_SMARTS_to_patterns(self.reactant1_dict)
        self.pattern_reactant2 = self.from_SMARTS_to_patterns(self.reactant2_dict)
        self.matched_reactants = None
        REACTANT1_PREFIX = self.fromReactionFullNameToReactantName("", 1)
        REACTANT2_PREFIX = self.fromReactionFullNameToReactantName("", 2)
        REACTIONS_NAMES = list(reaction_smarts.keys())

    def fromReactionFullNameToReactantName(self, reactionName, reactantNum):
        return ("react%d_" + reactionName.replace(" ", "_")) % reactantNum

    def fromReactantNameToReactionFullName(self, reactantName):
        return reactantName.replace("react1", "").replace("react2", "")

    def from_SMARTS_to_patterns(self, smarts_dict):
        return tuple([(key, Chem.MolFromSmarts(val)) for key, val in smarts_dict.items() if val is not None])

    def assign_reactants_w_rxn_smarts(self, reactant_attach_ids: Dict[Chem.Mol, List[Tuple[int, int]]],
                                      reaction_name: str) -> Dict[str, Tuple[Chem.Mol, List[int], str]]:
        """
        This function is used to assign the reactant number to input reactants using the reaction SMARTS. For now it only
        works for bimolecular reactions.

        :return:
        """
        if len(reactant_attach_ids) == 1: # Performing for one reactant
            patt = self.reactant1_dict[reaction_name]
            attach_ids = set(next(iter(reactant_attach_ids.values()))) # get first and only value
            reactant_mol = next(iter(reactant_attach_ids.keys())) # get first and only key
            self.matched_reactants: Dict[str, Tuple[Chem.Mol, List[int], str]] = (
                self.format_matched_reactant_for_one(reactant_mol, attach_ids, patt))
            assert len(self.matched_reactants) != 0 , "Reactant could not be matched to only SMARTS in reaction."
            return self.matched_reactants
        r1: Chem.Mol = list(reactant_attach_ids.keys())[0]
        r2: Chem.Mol = list(reactant_attach_ids.keys())[1]
        # Check to make sure they are not the same
        similarity = DataStructs.FingerprintSimilarity(rdMolDescriptors.GetMorganFingerprintAsBitVect(r1, 2),
                                                       rdMolDescriptors.GetMorganFingerprintAsBitVect(r2, 2))
        if similarity == 1.0:
            raise ValueError("The two reactants are the same.")
            return None
        r1_attach_ids: List[Tuple[int, int]] = set(reactant_attach_ids[r1])
        r2_attach_ids: List[Tuple[int, int]] = set(reactant_attach_ids[r2])
        patt1 = self.reactant1_dict[reaction_name]
        patt2 = self.reactant2_dict[reaction_name]
        self.matched_reactants = {patt1: None, patt2: None}
        found_1 = self.find_matching_atoms(r1, patt1, patt2, r1_attach_ids)
        found_2 = self.find_matching_atoms(r2, patt1, patt2, r2_attach_ids)
        # Check that both reactants have been found, and that they are not the same
        assert self.check_found_reactants(found_1, found_2, reaction_name, r1, r2)
        return self.matched_reactants

    def check_found_reactants(self, found_1: Dict[str, List], found_2: Dict[str, List], reaction_name: str,
                              r1: Chem.Mol, r2: Chem.Mol) -> bool:
        """
        This function checks that both reactants have been found, and that they are not the same.
        """
        if not found_1["r1"] and not found_2["r1"]:
            raise ValueError(f"The reactants do not match the reaction SMARTS in reaction {reaction_name} in "
                             f"mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}.")
            return None
        if found_1["r1"] and found_2["r1"]:
            raise ValueError(f"The reactants are the same in reaction {reaction_name} in "
                             f"mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}.")
            return None
        if not found_1["r2"] and not found_2["r2"]:
            print(
                f'WARNING: No atoms found involved in reaction {reaction_name} in '
                f'mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}')
            return None
        if found_1["r2"] and found_2["r2"]:
            print(
                f'WARNING: Both reactants have atoms involved in reaction {reaction_name} in '
                f'mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}')
            return None
        return True

    def find_matching_atoms(self, reactant: Chem.Mol, patt1: str, patt2: str, attachment_idx: set) -> (
            Dict)[str, Tuple[Chem.Mol, List[int], str]]:
        """
        This function finds the matched atoms in a reactant against both reactant SMARTS.
        """
        mol_patt1 = Chem.MolFromSmarts(patt1)
        mol_patt2 = Chem.MolFromSmarts(patt2)
        found = {'r1': False, 'r2': False}
        for mol_pattern, str_pattern, patt in zip([mol_patt1, mol_patt2], [patt1, patt2], ['r1', 'r2']):
            if found[patt]:
                continue
            matches = reactant.GetSubstructMatches(mol_pattern)
            for match in matches:
                if attachment_idx.intersection(match):
                    self.matched_reactants[str_pattern] = (reactant, match, patt)
                    found[patt] = match
                    break
        return found

    def format_matched_reactant_for_one(self, reactant_mol: Chem.Mol, attach_ids: List[Tuple[int, int]], patt: str) -> (
            Dict)[str, Tuple[Chem.Mol, List[int], str]]:
        """
        Formats matched reactants for one reactant.
        """
        matched_reactants = {}
        mol_patt = Chem.MolFromSmarts(patt)
        matches = reactant_mol.GetSubstructMatches(mol_patt)
        for match in matches:
            if attach_ids.intersection(match):
                matched_reactants[patt] = (reactant_mol, match, 'r1')
                break
            # Lax check if no attachment ids are found
            if len(match) == mol_patt.GetNumAtoms():
                matched_reactants[patt] = (reactant_mol, match, 'r1')
                break
        return matched_reactants



    

