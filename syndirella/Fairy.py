#!/usr/bin/env python3
"""
syndirella.cobblers_workshop.Fairy.py

This class 'Fairy' is used to provide reactants that offer additional similar cheaper reactants or filter out
reactants based on simple filters.
"""
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple, Optional)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import json
from syndirella.cli_defaults import cli_default_settings
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

class Fairy:
    def __init__(self):
        with open(cli_default_settings['reactant_filters_path']) as f:
            # the reactant_filters has structure of
            # {'filter_name': 'filter_smarts'}
            self.reactant_filters: Dict[str, str] = json.load(f)
        with open(cli_default_settings['additional_rxn_options_path']) as f:
            # the additional_rxn_options has structure of
            # {'reaction_name': 'additional_reaction_name'}
            self.additional_rxn_options: Dict[str, str] = json.load(f)

    def do_i_need_alterative_route(self, reaction_names: List[str]) -> bool:
        """
        This function is used to check if the reaction names need an alternative route.
        """
        for reaction_name in reaction_names:
            if reaction_name in self.additional_rxn_options:
                return True
        return False

    def find(self, reactant: Chem.Mol, reaction_name: str) -> List[str]:
        """
        This function is used to find additional similar reactants that are cheaper as defined in the
        self.reactant_filters.
        """
        orig_reactant: str = Chem.MolToSmiles(reactant)
        similar_reactant: str = None
        try:
            # Attempt to get the reactant filters using the reaction name
            reactant_filters: Dict[str, Any] = self.reactant_filters[reaction_name]
            print(f"Reaction name '{reaction_name}' found in reactant filters. Getting cheaper reactants...")
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
        # need to do validity check
        return similar_reactant

    def filter(self, hits: List[Tuple[str, float]]) -> Dict[str, float]:
        """
        Main entry to class to filter out reactants.
        """
        # Convert hits to mols
        mols = [Chem.MolFromSmiles(info[0]) for info in hits]
        print(f'Found {len(mols)} before filtering.')
        # Do simple filters
        mols = self.simple_filters(mols)
        # Print the difference
        self.print_diff(hits, mols, "simple filters")
        # Return dictionary
        return self._format_for_output(mols, hits)
        return {Chem.MolToSmiles(mol): hits[Chem.MolToSmiles(mol)] for mol in mols}

    def _format_for_output(self, mols: List[Chem.Mol], hits: List[Tuple[str, float]]) -> Dict[
        str, float]:
        """
        Formats the output of the filtered mols to a dictionary where the key is the SMILES of the molecule in the
        mols list, and the value is the float value that is paired with that SMILES in the hits tuple list.
        Comparison is done by Tanimoto similarity.

        Args:
        - mols: List of RDKit Mol objects.
        - hits: List of tuples, where each tuple contains a SMILES string and a float value.

        Returns:
        - A dictionary with SMILES strings as keys and float values as values.
        """
        hits_info: Dict[str, float] = {}
        # Precompute fingerprints for the hits for efficiency
        hits_fps = [(Chem.MolFromSmiles(smiles), score) for smiles, score in hits]
        hits_fps = [(AllChem.GetMorganFingerprintAsBitVect(mol, 2), score) for mol, score in hits_fps if
                    mol is not None]
        for mol in mols:
            mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            for hit_fp, value in hits_fps:
                similarity = TanimotoSimilarity(mol_fp, hit_fp)
                if similarity == 1.0:
                    mol_smiles = Chem.MolToSmiles(mol)
                    hits_info[mol_smiles] = value
                    break  # Assuming only one match is needed, break after the first match is found
        return hits_info

    def simple_filters(self, mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """
        This function is used to filter out mols that have
            chirality specification
            repeat mols
            non-abundant isotopes
        """
        # Remove mols with chirality specification, CHI_UNSPECIFIED is the non-chiral specification
        mols = [mol for mol in mols if all(atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED
                                           for atom in mol.GetAtoms())]
        # Remove repeat SMILES
        mols = Fairy.remove_repeat_mols(mols)
        # Remove non-abundant isotopes
        mols = Fairy._remove_non_abundant_isotopes(mols)
        return mols

    @staticmethod
    def remove_repeat_mols(mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """
        Converts mols to fingerprints, checks for Tanimoto similarity of 1, and removes repeats.
        """
        print("Removing repeat analogues...")
        unique_mols = []
        fingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048) for mol in mols]
        seen = set()  # To keep track of indices of molecules that are duplicates
        for i, fp1 in enumerate(fingerprints):
            if i not in seen:
                unique_mols.append(mols[i])
                for j, fp2 in enumerate(fingerprints[i + 1:], start=i + 1):
                    similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
                    if similarity == 1:
                        seen.add(j)
        return unique_mols

    @staticmethod
    def _remove_non_abundant_isotopes(mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """
        This function is used to remove mols that have non-abundant isotopes.
        """
        mols = [mol for mol in mols if
                not any((atom.GetIsotope() not in [0, 12] and atom.GetAtomicNum() == 6) or # Carbon
                        (atom.GetIsotope() not in [0, 14] and atom.GetAtomicNum() == 7) or # Nitrogen
                        (atom.GetIsotope() not in [0, 1] and atom.GetAtomicNum() == 1) or # Hydrogen
                        (atom.GetIsotope() not in [0, 16] and atom.GetAtomicNum() == 8) # Oxygen
                        for atom in mol.GetAtoms())]
        return mols

    @staticmethod
    def _remove_chirality_from_mol(mol: Chem.Mol) -> Chem.Mol:
        """
        This function is used to remove chirality from a single molecule.
        """
        [atom.SetChiralTag(Chem.CHI_UNSPECIFIED) for atom in mol.GetAtoms()]
        return mol

    @staticmethod
    def remove_chirality(mols: List[Chem.Mol]) -> List[Chem.Mol]:
        """
        This function is used to remove mols that have chirality specification.
        """
        print("Removing chirality from analogues...")
        mols = [Fairy._remove_chirality_from_mol(mol) for mol in mols]
        return mols


    def print_diff(self, mols: List[Chem.Mol], valid_mols: List[Chem.Mol], desc: str) -> None:
        """
        This function is used to print the difference between the original number of analogues and the number of
        valid analogues.
        """
        assert len(valid_mols) <= len(mols), ("Problem with finding valid molecules. There are more than were in "
                                              "the original list of molecules.")
        num_filtered = len(mols) - len(valid_mols)
        percent_diff = round((num_filtered / len(mols)) * 100, 2)
        print(f'Removed {num_filtered} molecules ({percent_diff}%) by {desc}.')

    def get_final_routes(self, routes: List[List[Dict]]) -> List[List[Dict]]:
        """
        This function returns the final routes to elaborate on. It is called by Cobbler class.
        """
        # check what reactions are in it
        reaction_names = [[reaction['name'].replace(" ","_") for reaction in route] for route in routes]
        # get first route
        first_route_names = reaction_names[0]
        if len(first_route_names) == 1:
            print(f"The route found is {len(first_route_names)} step. The forward synthesis is:")
            print(f"{first_route_names[::-1]}") # reverse the order to make it forward synthesis
        else:
            print(f"The first route found is {len(first_route_names)} steps. The forward synthesis is:")
            print(f"{first_route_names[::-1]}") # reverse the order to make it forward synthesis
        final_routes: List[List[Dict]] = self.get_additional_routes(first_route_names, routes)
        return final_routes

    def get_additional_routes(self, first_route_names: List[str], routes: List[List[Dict]]) -> List[List[Dict]]:
        """
        This function is used to get additional routes if specified in fairy filters.
        """
        additional_route = None
        # look at each reaction in the first route
        for reaction_name in first_route_names:
            # get additional reactions if specified in fairy filters
            if reaction_name in self.additional_rxn_options:
                additional_reaction_name = self.additional_rxn_options[reaction_name]
                print(f"Additional reaction for '{reaction_name}' found in fairy filters. "
                      f"Getting additional routes containing '{additional_reaction_name}'...")
                # just get the first additional route
                additional_route: List[List[Dict]] = self.get_additional_route(additional_reaction_name, routes)
        final_route: List[List[Dict]] = [routes[0]]
        assert final_route != additional_route, "The first route and the additional route are the same."
        if additional_route is not None:
            final_route.append(additional_route)
        return final_route

    def get_additional_route(self, additional_reaction_name: str, routes: List[List[Dict]]) -> List[List[Dict]]:
        """
        This function is used to get additional routes containing the additional reaction name.
        """
        # get additional route that is not the first one containing the additional reaction name
        additional_routes = [route for route in routes if
                             any(additional_reaction_name in reaction['name'] for reaction in route)]
        return additional_routes[0] # return the first additional route

