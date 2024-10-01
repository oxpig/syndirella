#!/usr/bin/env python3
"""
syndirella.route.fairy.py

This module provides functions to find similar cheaper reactants, filter out
reactants based on simple filters, fingerprint generation, and others...
"""

from typing import Any, List, Dict, Tuple, Optional
from rdkit.Chem import AllChem, rdinchi
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import json
from syndirella.cli_defaults import cli_default_settings
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
import logging
from rdkit.Chem import rdFingerprintGenerator
from syndirella.error import *

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
    """Filter out molecules based on repeats, and non-abundant isotopes."""
    mols = remove_repeat_mols(mols)
    mols = remove_non_abundant_isotopes(mols)
    mols = remove_hydrogen_ions(mols)
    return mols

def remove_hydrogen_ions(mols: List[Chem.Mol]) -> List[Chem.Mol]:
    """Remove molecules with only hydrogen ions."""
    for i, mol in enumerate(mols):
        mol = Chem.RWMol(mol) # make editable
        atoms_to_remove = []
        for atom in mol.GetAtoms():
            # remove single hydrogen ions
            if atom.GetAtomicNum() == 1 and atom.GetHybridization() == Chem.HybridizationType.S:
                atoms_to_remove.append(atom.GetIdx())
        try:
            for atom_idx in atoms_to_remove[::-1]: # have to reverse it for some reason
                mol.RemoveAtom(atom_idx)
        except Exception as e: # don't know exact Error here
            logger.info(f"Error removing hydrogen ion: {e}")
            continue
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

def calculate_tanimoto(mol1: Chem.Mol,
                       mol2: Chem.Mol) -> float:
    if mol1 is None or mol2 is None:
        return 0
    # do full circle back to smiles to mol
    mol1, mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1)), Chem.MolFromSmiles(Chem.MolToSmiles(mol2))
    fp1 = get_morgan_fingerprint(mol1)
    fp2 = get_morgan_fingerprint(mol2)
    return TanimotoSimilarity(fp1, fp2)


def find_reaction_by_name(reaction_name: str, additional_rxn_options: List[Dict[str, str]]) -> Optional[Dict[str, Any]]:
    """
    Find a reaction by name in the additional_rxn_options.
    """
    for reaction in additional_rxn_options:
        if reaction['name'] == reaction_name:
            return reaction
    return None

def get_morgan_fingerprint(mol: Chem.Mol) -> rdFingerprintGenerator.MorganFP:
    """
    Get the fingerprint of a molecule.
    """
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return mfpgen.GetFingerprint(mol)


def generate_inchi_ID(smiles: str | None = None, mol: Chem.Mol | None = None, isomeric: bool = False) -> str | MolError:
    """
    This function is used to generate a unique id for the route just using the scaffold.
    """
    if smiles is not None and isomeric is False:
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=False) # remove isomeric information
        mol = Chem.MolFromSmiles(smiles)
    if smiles is not None and Chem.MolFromSmiles(smiles) is None:
        logger.critical(f"Could not create a molecule from the smiles {smiles}.")
        return MolError(message=f"Could not create a molecule from the smiles {smiles}.",
                        smiles=smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles)
    ID = rdinchi.MolToInchi(mol)
    id = rdinchi.InchiToInchiKey(ID[0])
    return id

def check_mol_sanity(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Check if the molecule can be sanitized.
    """
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return mol
    except ValueError:
        return None

