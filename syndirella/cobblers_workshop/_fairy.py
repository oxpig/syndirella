#!/usr/bin/env python3
"""
syndirella.cobblers_workshop._fairy.py

This class 'Fairy' is used to provide reactants that offer additional similar cheaper reactants or filter out
reactants based on simple filters.
"""
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple, Optional)
from rdkit import Chem
import json
from syndirella._cli_defaults import cli_default_settings
from syndirella.cobblers_workshop._reaction import Reaction

class Fairy:
    def __init__(self, reaction_name: str):
        with open(cli_default_settings['reactant_filters_path']) as f:
            # the reactant_filters has structure of
            # {'filter_name': 'filter_smarts'}
            self.reactant_filters: Dict[str, str] = json.load(f)
        self.reaction_name = reaction_name

    def find(self, reactant: Chem.Mol) -> List[str]:
        """
        This function is used to find additional similar reactants that are cheaper as defined in the
        self.reactant_filters.
        """
        orig_reactant: str = Chem.MolToSmiles(reactant)
        try:
            # Attempt to get the reactant filters using the reaction name
            reactant_filters: Dict[str, Any] = self.reactant_filters[self.reaction_name]
            print(f"Reaction name '{self.reaction_name}' found in reactant filters. Getting cheaper reactants...")
            # get similar reactant
            similar_reactant: str = self._get_similar_reactant(reactant, reactant_filters)
        except KeyError:
            # Handle the case where the reaction name is not a key in the reactant_filters
            reactant_filters = None
        if similar_reactant is None or reactant_filters is None: # if no similar reactant is found
            return [orig_reactant]
        return [orig_reactant, similar_reactant]

    def _get_similar_reactant(self, reactant: Chem.Mol, reactant_filters: Dict[str, str]) -> str:
        """
        Makes similar reactant to query reactant from the reactant_filters.
        """
        # Initialize a variable to store the details if exactly one match is found
        matched_details = None
        # smarts match reactant to filter
        for filter_smarts, to_add_details in reactant_filters.items():
            filter_smarts: str
            to_add_details: Dict[str, Any]
            filter_mol = Chem.MolFromSmarts(filter_smarts)
            # make sure reactant has substruct match to only 1 filter
            if reactant.HasSubstructMatch(filter_mol):
                if matched_details is None:
                    # This is the first match, store the details
                    to_replace = filter_smarts
                    matched_details = to_add_details
                else:
                    # Found another match, meaning more than one match is found
                    # Reset matched_details to indicate multiple matches
                    matched_details = None
                    break  # Exit the loop as we found more than one match
        # Proceed only if exactly one match was found
        if matched_details is not None:
            # make similar reactant with the details that matched exactly one filter
            similar_reactant = self._make_similar_reactant(reactant, to_replace, matched_details)
            return similar_reactant
        else:
            # If no match or more than one match is found, return None
            return None

    def _make_similar_reactant(self, reactant: Chem.Mol, to_replace: str, matched_details: Dict[str, Any]) -> str:
        """
        Makes similar reactant from SMARTS to add and replace part of original reactant.
        """
        name = matched_details['name']
        print(f"Performing '{name}'...")
        # Get the details to add
        to_add: str = matched_details['to_add']
        connecting_atom_id: int = matched_details['connecting_atom_id']
        to_replace = Chem.MolFromSmarts(to_replace)
        to_add = Chem.MolFromSmarts(to_add)
        # Replace the part of the reactant with to_replace with to_add
        similar_reactant = Chem.ReplaceSubstructs(reactant, to_replace, to_add,
                                                  replacementConnectionPoint=connecting_atom_id, replaceAll=True)
        # TODO: Could add functionality to check that the similar reactant was changed at correct point
        # Convert the similar reactant to smiles
        similar_reactant = Chem.MolToSmiles(similar_reactant[0])
        return similar_reactant

    def filter(self, hits: Dict[str, float]) -> Dict[str, float]:
        """
        Main entry to class to filter out reactants.
        """
        # Convert hits to mols
        mols = [Chem.MolFromSmiles(smile) for smile in hits.keys()]
        # Do simple filters
        mols = self.simple_filters(mols)
        # Return dictionary
        return {Chem.MolToSmiles(mol): hits[Chem.MolToSmiles(mol)] for mol in mols}

    def simple_filters(self, mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """
        This function is used to filter out mols that have
            chirality specification
            repeat mols
            non-abundant isotopes
        """
        # Remove mols with chirality specification
        mols = [mol for mol in mols if not mol.HasProp('_ChiralityPossible')]
        # Remove repeat mols
        mols = list(set(mols))
        # Remove non-abundant isotopes
        mols = [mol for mol in mols if not any(atom.GetIsotope() not in [0, 12, 13, 14, 15] for atom in mol.GetAtoms())]
        return mols

    def print_diff(self, mols: List[Chem.Mol], valid_mols: List[Chem.Mol], name: str) -> None:
        """
        This function is used to print the difference between the original number of analogues and the number of
        valid analogues.
        """
        assert len(valid_mols) <= len(mols), ("Problem with finding valid molecules. There are more than were in "
                                              "the original list of molecules.")
        num_filtered = len(mols) - len(valid_mols)
        percent_diff = round((num_filtered / len(mols)) * 100, 2)
        print(f'Removed {num_filtered} invalid molecules ({percent_diff}%) that did not have {name} substructure.')

