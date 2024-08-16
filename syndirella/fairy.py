#!/usr/bin/env python3
"""
syndirella.cobblers_workshop.fairy.py

This module provides functions to find similar cheaper reactants or filter out
reactants based on simple filters.
"""

from typing import Any, List, Dict, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import json
from syndirella.cli_defaults import cli_default_settings
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
import logging
from rdkit.Chem import rdFingerprintGenerator

# Set up logger
logger = logging.getLogger(__name__)


# Load necessary data
def load_reactant_filters() -> Dict[str, Dict[str, Any]]:
    with open(cli_default_settings['reactant_filters_path']) as f:
        return json.load(f)


def load_additional_rxn_options() -> List[Dict[str, str]]:
    with open(cli_default_settings['additional_rxn_options_path']) as f:
        return json.load(f)


def do_i_need_alternative_route(reaction_names: List[str], additional_rxn_options: List[Dict[str, str]]) -> bool:
    """
    Check if the reaction names need an alternative route.
    """
    return any(reaction_name in additional_rxn_options for reaction_name in reaction_names)


def find_similar_reactants(reactant: Chem.Mol, reaction_name: str, reactant_filters: Dict[str, Dict[str, Any]]) -> List[
    str] | None:
    """
    Find additional similar reactants that are cheaper as defined in the reactant_filters.
    """
    orig_reactant: str = Chem.MolToSmiles(reactant)
    similar_reactant: Optional[str] = None
    try:
        filters_for_reaction: Dict[str, Any] = reactant_filters[reaction_name]
        logger.info(f"Reaction name '{reaction_name}' found in reactant filters. Getting cheaper reactants...")
        similar_reactant = get_similar_reactant(reactant, filters_for_reaction)
    except KeyError:
        filters_for_reaction = None

    if similar_reactant is None or filters_for_reaction is None:
        return [orig_reactant]

    return [orig_reactant, similar_reactant]


def get_similar_reactant(reactant: Chem.Mol, reactant_filters: Dict[str, str]) -> Optional[str]:
    """Make similar reactant to query reactant from the reactant_filters."""
    matched_details: Optional[Dict[str, Any]] = None

    for filter_smarts, to_add_details in reactant_filters.items():
        filter_mol = Chem.MolFromSmarts(filter_smarts)
        if reactant.HasSubstructMatch(filter_mol):
            if matched_details is None:
                to_replace = filter_smarts
                matched_details = to_add_details
            else:
                matched_details = None
                break

    if matched_details is not None:
        similar_reactant = make_similar_reactant(reactant, to_replace, matched_details)
        return similar_reactant
    return None


def make_similar_reactant(reactant: Chem.Mol, to_replace: str, matched_details: Dict[str, Any]) -> str:
    """Create a similar reactant by modifying the original reactant according to matched_details."""
    name = matched_details['name']
    logger.info(f"Performing '{name}'...")

    to_add = Chem.MolFromSmarts(matched_details['to_add'])
    connecting_atom_id = matched_details['connecting_atom_id']
    to_replace_mol = Chem.MolFromSmarts(to_replace)

    similar_reactant = Chem.ReplaceSubstructs(reactant, to_replace_mol, to_add,
                                              replacementConnectionPoint=connecting_atom_id, replaceAll=True)
    similar_reactant_smiles = Chem.MolToSmiles(similar_reactant[0])

    return similar_reactant_smiles


def filter_molecules(hits: List[Tuple[str, float]]) -> Dict[str, float]:
    """
    Filter out reactants and return a dictionary with SMILES strings as keys and their associated scores as values.
    """
    mols = [Chem.MolFromSmiles(info[0]) for info in hits]
    logger.info(f'Found {len(mols)} before filtering.')

    mols = simple_filters(mols)
    print_diff(hits, mols, "simple filters")

    return format_for_output(mols, hits)


def format_for_output(mols: List[Chem.Mol], hits: List[Tuple[str, float]]) -> Dict[str, float]:
    """
    Format the output of the filtered molecules to a dictionary with SMILES strings as keys and their scores as values.
    """
    hits_info: Dict[str, float] = {}
    hits_fps = [(Chem.MolFromSmiles(smiles), score) for smiles, score in hits]
    hits_fps = [(AllChem.GetMorganFingerprintAsBitVect(mol, 2), score) for mol, score in hits_fps if mol is not None]

    for mol in mols:
        mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        for hit_fp, value in hits_fps:
            similarity = TanimotoSimilarity(mol_fp, hit_fp)
            if similarity == 1.0:
                mol_smiles = Chem.MolToSmiles(mol)
                hits_info[mol_smiles] = value
                break

    return hits_info


def simple_filters(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    """Filter out molecules based on chirality, repeats, and non-abundant isotopes."""
    mols = [mol for mol in mols if all(atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED
                                       for atom in mol.GetAtoms())]
    mols = remove_repeat_mols(mols)
    mols = remove_non_abundant_isotopes(mols)
    return mols


def remove_repeat_mols(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    """
    Remove repeated molecules by checking for Tanimoto similarity of 1.
    """
    logger.info("Removing repeat analogues...")
    unique_mols = []
    fingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048) for mol in mols]
    seen = set()

    for i, fp1 in enumerate(fingerprints):
        if i not in seen:
            unique_mols.append(mols[i])
            for j, fp2 in enumerate(fingerprints[i + 1:], start=i + 1):
                similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
                if similarity == 1:
                    seen.add(j)

    return unique_mols


def remove_non_abundant_isotopes(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    """Remove molecules with non-abundant isotopes."""
    return [mol for mol in mols if not any(
        (atom.GetIsotope() not in [0, 12] and atom.GetAtomicNum() == 6) or
        (atom.GetIsotope() not in [0, 14] and atom.GetAtomicNum() == 7) or
        (atom.GetIsotope() not in [0, 1] and atom.GetAtomicNum() == 1) or
        (atom.GetIsotope() not in [0, 16] and atom.GetAtomicNum() == 8)
        for atom in mol.GetAtoms())]


def remove_chirality(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    """
    Remove chirality from all molecules in the list.
    """
    logger.info("Removing chirality from analogues...")
    return [remove_chirality_from_mol(mol) for mol in mols]


def remove_chirality_from_mol(mol: Chem.Mol) -> Chem.Mol:
    """Remove chirality from a single molecule."""
    for atom in mol.GetAtoms():
        atom.SetChiralTag(Chem.CHI_UNSPECIFIED)
    return mol


def print_diff(hits: List[Tuple[str, float]], valid_mols: List[Chem.Mol], desc: str) -> None:
    """Print the difference between the original number of molecules and the number of valid molecules."""
    num_filtered = len(hits) - len(valid_mols)
    percent_diff = round((num_filtered / len(hits)) * 100, 2)
    logger.info(f'Removed {num_filtered} molecules ({percent_diff}%) by {desc}.')


def get_final_routes(routes: List[List[Dict]], additional_rxn_options: Dict[str, str]) -> List[List[Dict]]:
    """
    Return the final routes to elaborate on.
    """
    reaction_names = [[reaction['name'].replace(" ", "_") for reaction in route] for route in routes]
    first_route_names = reaction_names[0]

    if len(first_route_names) == 1:
        logger.info(f"The route found is {len(first_route_names)} step. "
                    f"The forward synthesis is: {first_route_names[::-1]}")
    else:
        logger.info(f"The first route found is {len(first_route_names)} steps. "
                    f"The forward synthesis is: {first_route_names[::-1]}")

    final_routes = get_additional_routes(first_route_names, routes, additional_rxn_options)
    return final_routes


def get_additional_routes(first_route_names: List[str], routes: List[List[Dict]],
                          additional_rxn_options: Dict[str, str]) -> List[List[Dict]]:
    """Get additional routes if specified in filters."""
    additional_route = None

    for reaction_name in first_route_names:
        if reaction_name in additional_rxn_options:
            additional_reaction_name = additional_rxn_options[reaction_name]
            logger.info(f"Additional reaction for '{reaction_name}' found in filters. "
                        f"Getting additional routes containing '{additional_reaction_name}'...")
            additional_route = [route for route in routes if
                                any(additional_reaction_name in reaction['name'] for reaction in route)]
            if additional_route:
                additional_route = additional_route[0]

    final_route = [routes[0]]
    if additional_route and additional_route != final_route:
        final_route.append(additional_route)

    return final_route

def calculate_tanimoto(mol1: Chem.Mol,
                       mol2: Chem.Mol) -> float:
    if mol1 is None or mol2 is None:
        return 0
    # do full circle back to smiles to mol
    mol1, mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1)), Chem.MolFromSmiles(Chem.MolToSmiles(mol2))
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp1 = mfpgen.GetFingerprint(mol1)
    fp2 = mfpgen.GetFingerprint(mol2)
    return TanimotoSimilarity(fp1, fp2)


# def define_reactions(workshop: CobblersWorkshop) -> CobblersWorkshop:
#     """
#     Define the reactions in the workshop.
#     """
#     for step in range(workshop.num_steps):
#         cobbler_bench = workshop.get_cobbler_bench(step)
#         cobbler_bench.define_reaction()
#     return workshop


def find_reaction_by_name(reaction_name: str, additional_rxn_options: List[Dict[str, str]]) -> Optional[Dict[str, Any]]:
    """
    Find a reaction by name in the additional_rxn_options.
    """
    for reaction in additional_rxn_options:
        if reaction['name'] == reaction_name:
            return reaction
    return None


def get_additional_reactants(bench, additional_reaction: Dict[str, Any]) -> Tuple[str, str]:
    """Edit the reactants."""
    # Implementation goes here


# def get_additional_route_from_workshop(workshop: CobblersWorkshop, additional_rxn_options: List[Dict[str, str]]) -> Optional[Any]:
#     """
#     Get an additional route from a workshop object where the reactions have all been defined and checked.
#     """
#     additional_route = None
#     reactants: List[Tuple[str, str]] = []
#     reaction_names: List[str] = []
#
#     for bench in workshop.cobbler_benches:
#         if bench.reaction_name in additional_rxn_options:
#             additional_reaction = find_reaction_by_name(bench.reaction_name, additional_rxn_options)
#             logger.info(f"Additional reaction for '{bench.reaction_name}' found in filters. "
#                         f"Getting additional routes containing '{additional_reaction['name']}'...")
#             new_reactants, new_reaction_name = get_additional_reactants(bench, additional_reaction)
#
#     return additional_route


# def get_additional_route(workshop: CobblersWorkshop, additional_rxn_options: List[Dict[str, str]]) -> CobblersWorkshop:
#     """
#     Get an alternative route if it contains specified reactions.
#     """
#     if do_i_need_alternative_route(workshop.reaction_names, additional_rxn_options):
#         logger.info(f"Found the need for an alternative route for compound {workshop.product}.")
#         workshop: CobblersWorkshop = define_reactions(workshop)
#         additional_route: CobblersWorkshop = get_additional_route_from_workshop(workshop, additional_rxn_options)
#         return additional_route
#     return None