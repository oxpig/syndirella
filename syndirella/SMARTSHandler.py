#!/usr/bin/env python3
"""
SMARTSHandler.py

This module contains the SMARTSHandler class. This class contains information about the reaction SMARTS.
"""
import json
import logging
from collections import OrderedDict
from rdkit import Chem, DataStructs
from rdkit.Chem import Mol
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from typing import (Set, List, Dict, Tuple, Any)
from .cli_defaults import cli_default_settings
import syndirella.fairy as fairy
from .error import SMARTSError


class SMARTSHandler:

    def __init__(self):
        with open(cli_default_settings['rxn_smarts_path']) as f:
            reaction_smarts = json.load(f)
        self.reaction_smarts = {name: ReactionFromSmarts(val) for name, val in reaction_smarts.items()}
        self.reactant1_dict = OrderedDict()
        self.reactant2_dict = OrderedDict()
        self.product_smarts = OrderedDict()
        self.n_reactants_per_reaction = OrderedDict()
        for name, smart in reaction_smarts.items():
            reactants, prod = smart.split(">>")
            try:
                react1, react2 = reactants.split(".")
            except ValueError:
                react1, react2 = reactants, None
            self.reactant1_dict[name] = react1
            self.reactant2_dict[name] = react2
            self.product_smarts[name] = prod
            self.n_reactants_per_reaction = {name: 1 if react2 is None else 2}
        self.pattern_products = self.from_SMARTS_to_patterns(self.product_smarts)
        self.pattern_reactant1 = self.from_SMARTS_to_patterns(self.reactant1_dict)
        self.pattern_reactant2 = self.from_SMARTS_to_patterns(self.reactant2_dict)
        self.matched_reactants = None
        self.found_1 = None
        self.found_2 = None
        REACTANT1_PREFIX = self.fromReactionFullNameToReactantName("", 1)
        REACTANT2_PREFIX = self.fromReactionFullNameToReactantName("", 2)
        REACTIONS_NAMES = list(reaction_smarts.keys())
        self.logger = logging.getLogger(f"{__name__}")

    def fromReactionFullNameToReactantName(self, reactionName, reactantNum):
        return ("react%d_" + reactionName.replace(" ", "_")) % reactantNum

    def fromReactantNameToReactionFullName(self, reactantName):
        return reactantName.replace("react1", "").replace("react2", "")

    def from_SMARTS_to_patterns(self, smarts_dict):
        return tuple([(key, Chem.MolFromSmarts(val)) for key, val in smarts_dict.items() if val is not None])

    def assign_reactants_w_rxn_smarts(self,
                                      product: Chem.Mol,
                                      reactant_attach_ids: Dict[Chem.Mol, List[Tuple[int, int]]],
                                      reaction_name: str) -> dict[
                                                                 str, tuple[Mol, list[int], str]] | None | dict[
                                                                 Any | None, None]:
        """
        This function is used to assign the reactant number to input reactants using the reaction SMARTS. For now it
        only supports bimolecular reactions.
        """
        if len(reactant_attach_ids) == 1:  # Performing for one reactant
            patt = self.reactant1_dict[reaction_name]
            attach_ids = set(next(iter(reactant_attach_ids.values())))  # get first and only value
            reactant_mol = next(iter(reactant_attach_ids.keys()))  # get first and only key
            self.matched_reactants: Dict[str, Tuple[Chem.Mol, List[int], str]] = (
                self.format_matched_reactant_for_one(reactant_mol, attach_ids, patt))
            if len(self.matched_reactants) == 0:
                message = "Reactant could not be matched to only reactant SMARTS in reaction."
                self.logger.critical(message)
                raise SMARTSError(mol=product, message=message)
            return self.matched_reactants
        r1: Chem.Mol = list(reactant_attach_ids.keys())[0]
        r2: Chem.Mol = list(reactant_attach_ids.keys())[1]
        # Check to make sure they are not the same
        similarity = DataStructs.FingerprintSimilarity(fairy.get_morgan_fingerprint(r1),
                                                       fairy.get_morgan_fingerprint(r2))
        if similarity == 1.0:
            message = "The two reactants are the same."
            self.logger.critical(message)
            raise SMARTSError(mol=product, message=message)
        r1_attach_ids: Set[Tuple[int, int]] = set(reactant_attach_ids[r1])
        r2_attach_ids: Set[Tuple[int, int]] = set(reactant_attach_ids[r2])
        patt1 = self.reactant1_dict[reaction_name]
        patt2 = self.reactant2_dict[reaction_name]
        self.matched_reactants = {patt1: None, patt2: None}
        found_1: Dict[str, Tuple[int] | bool] = self.find_matching_atoms(r1, patt1, patt2, r1_attach_ids)
        found_2: Dict[str, Tuple[int] | bool] = self.find_matching_atoms(r2, patt1, patt2, r2_attach_ids)
        # Check that both reactants have been found
        if not self.check_found_reactants(product, found_1, found_2, reaction_name, r1, r2):
            # will return False when both reactants match both reactant SMARTS. i.e. formation of urea.
            self.seperate_matching_reactants(r1, r2, found_1, found_2, patt1, patt2)
        self.found_1 = found_1
        self.found_2 = found_2
        return self.matched_reactants

    def seperate_matching_reactants(self,
                                    r1: Chem.Mol,
                                    r2: Chem.Mol,
                                    found_1: Dict[str, bool],
                                    found_2: Dict[str, bool],
                                    patt1: str,
                                    patt2: str) -> Dict[str, Tuple[Chem.Mol, List[int], str]]:
        """
        This function is used to fix edge cases:
        - to separate reactants that match both reactant SMARTS.
        - one reactant matches both reactant SMARTS.

        These found_1 and found_2 dictionaries are a bit confusing:
        Key: 'r1' or 'r2'
        Value: False if no match, List of atom indices within reactant that match the SMARTS of that reactant.
        """
        # Both reactants match both reactant SMARTS
        if found_1["r1"] and found_2["r1"] and found_2["r2"] and found_1["r2"]:
            self.matched_reactants[patt1] = (r1, found_1["r1"], 'r1')
            self.matched_reactants[patt2] = (r2, found_2["r2"], 'r2') # make found_2 r2
        # reactant 1 matches both reactant SMARTS
        elif found_1["r1"] and found_1["r2"]:
            if found_2["r2"]: # change r2 to found_2
                self.matched_reactants[patt1] = (r1, found_1["r1"], 'r1')
                self.matched_reactants[patt2] = (r2, found_2["r2"], 'r2')
            if found_2["r1"]: # change r1 to found_2
                self.matched_reactants[patt1] = (r2, found_2["r1"], 'r1')
                self.matched_reactants[patt2] = (r1, found_1["r2"], 'r2')
        # reactant 2 matches both reactant SMARTS
        elif found_2["r1"] and found_2["r2"]:
            if found_1["r2"]: # change r2 to found_1
                self.matched_reactants[patt1] = (r2, found_2["r1"], 'r1')
                self.matched_reactants[patt2] = (r1, found_1["r2"], 'r2')
            if found_1["r1"]: # change r1 to found_1
                self.matched_reactants[patt1] = (r1, found_1["r1"], 'r1')
                self.matched_reactants[patt2] = (r2, found_2["r1"], 'r2')

    def check_found_reactants(self,
                              product: Chem.Mol,
                              found_1: Dict[str, bool],
                              found_2: Dict[str, bool],
                              reaction_name: str,
                              r1: Chem.Mol,
                              r2: Chem.Mol) -> bool:
        """
        This function checks that both reactants have been found. It raises a warning if both reactants are matched
        to the same reactant SMARTS.

        Returns False when both reactants match both reactant SMARTS.
        """
        if not found_1["r1"] and not found_2["r1"]:
            self.logger.critical(f"The reactants do not match the reaction SMARTS in reaction {reaction_name} in "
                                 f"mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}.")
            raise SMARTSError(message=f"The reactants do not match the reaction SMARTS in reaction {reaction_name} in "
                                      f"mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}.", mol=product)

        if not found_1["r2"] and not found_2["r2"]:
            self.logger.critical(
                f"No atoms found involved in reaction {reaction_name} in "
                f"mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}")
            raise SMARTSError(message=
                              f"No atoms found involved in reaction {reaction_name} in "
                              f"mol {Chem.MolToSmiles(r1)} and {Chem.MolToSmiles(r2)}",
                              mol=product)

        if found_1["r1"] and found_2["r1"]:
            self.logger.warning(
                f"Both reactants ({Chem.MolToSmiles(r1)} {Chem.MolToSmiles(r2)}) have atoms found in SMARTS of 1st "
                f"reactant for {reaction_name}. This might cause selectivity issues downstream, but continuing.")
            if found_1["r2"] and found_2["r2"]:
                self.logger.warning(
                    f"Both reactants ({Chem.MolToSmiles(r1)} {Chem.MolToSmiles(r2)}) have atoms found in SMARTS of 2nd "
                    f"reactant for {reaction_name}. This might cause selectivity issues downstream, but continuing.")
            return False

        return True

    def find_matching_atoms(self, reactant: Chem.Mol, patt1: str, patt2: str, attachment_idx: set) -> Dict[str, Tuple[int] | bool]:
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
