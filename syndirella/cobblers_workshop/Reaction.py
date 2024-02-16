#!/usr/bin/env python3
"""
syndirella.cobblers_workshop.Reaction.py

This module contains the Reaction class. One instance of this object is used to describe a single reaction.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS
from typing import (List, Dict, Tuple)
from syndirella.constants.constants import DUMMY_SYMBOL
from syndirella.SMARTSHandler import SMARTSHandler
from syndirella.error import ReactionError

class Reaction():
    """
    This class contains information about the reaction. It is used to find the atoms on the reactants involved in
    the reaction, check the reaction atoms are the ones that are supposed to be connected.
    """
    def __init__(self, product: Chem.Mol, reactants: List[Chem.Mol], reaction_name: str, smarts_handler: SMARTSHandler):
        self.product: Chem.Mol = product
        self.reactants: List[Chem.Mol] = reactants
        self.reaction_name: str = reaction_name
        self.smarts_handler: SMARTSHandler = smarts_handler
        self.reaction_pattern: Chem.rdChemReactions = self.smarts_handler.reaction_smarts[self.reaction_name]
        self.all_attach_ids: List[int] = None
        self.matched_smarts_to_reactant: Dict[str, Tuple[Chem.Mol, List[int]]] = None

    def _replace_carboxylic_acid_hydroxy_with_dummy(self, mol: Chem.Mol) -> Chem.Mol:
        # Get the oxygen atoms in the molecule
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 8:  # Atomic number for Oxygen
                # Check if the oxygen is in a carboxylic acid group
                if len(atom.GetNeighbors()) == 1 and atom.GetNeighbors()[0].GetSymbol() == 'C':
                    carboxylic_carbon = atom.GetNeighbors()[0]
                    if len([nbr for nbr in carboxylic_carbon.GetNeighbors() if nbr.GetSymbol() == 'O']) == 2:
                        atom.SetAtomicNum(0)  # Set atomic number to 0, which is a dummy atom in RDKit
        return mol

    def _replace_halide_with_dummy(self, mol: Chem.Mol, atom_to_check_by_symbol=None) -> Chem.Mol:
        # Convert the symbol to atomic number
        if atom_to_check_by_symbol is not None:
            atom_to_check_num = Chem.GetPeriodicTable().GetAtomicNumber(atom_to_check_by_symbol)
        # Start editing the molecule
        rw_mol = Chem.RWMol(mol)
        for atom in rw_mol.GetAtoms():
            if atom.GetAtomicNum() in [9, 17, 35, 53]:  # For F, Cl, Br, I
                for neighbor in atom.GetNeighbors():
                    if atom_to_check_by_symbol is not None:
                        if neighbor.GetNum() == atom_to_check_num:
                            atom.SetAtomicNum(0)
                    else:
                        atom.SetAtomicNum(0)  # Set atomic number to 0 for dummy atom
                        break
        # Convert RWMol back to Mol object
        return rw_mol.GetMol()

    def _find_attachment_id_from_dummy(self, reactant: Chem.Mol, dummy_symbol="*") -> Tuple[List[int], List[List[int]]]:
        neig_idx_list = []
        dummy_idx_list = []
        for idx in range(reactant.GetNumAtoms()):
            if reactant.GetAtomWithIdx(idx).GetSymbol() == dummy_symbol:
                dummy_idx = idx
                neig_idx = [atom.GetIdx() for atom in reactant.GetAtomWithIdx(idx).GetNeighbors()]
                dummy_idx_list.append(dummy_idx)
                neig_idx_list.append(neig_idx)

        return dummy_idx_list, neig_idx_list

    def use_fmcs(self, reactant: Chem.Mol) -> Dict[int, int]:
        mcs = rdFMCS.FindMCS([reactant, self.product])
        mcs_smarts = Chem.MolFromSmarts(mcs.smartsString)
        product_matches = self.product.GetSubstructMatches(mcs_smarts)
        # Can return multiple product_matches, only take the first one
        product_match = product_matches[0]
        reactant_matches = reactant.GetSubstructMatches(mcs_smarts)
        reactant_match = reactant_matches[0]
        product_to_reactant_mapping = dict(zip(product_match, reactant_match))
        return product_to_reactant_mapping

    def find_attachment_id_for_reactant(self, reactant: Chem.Mol) -> List[Tuple[int, int]]:
        """
        This function is used to find the attachment indices of a single reactant in the reaction.
        :param reactant: a single reactant molecule
        :returns a list of tuples (attachmentIdx_whole, attachmentIdx_subMol)
        """
        # product_to_reactant_mapping = self.product.GetSubstructMatch(reactant)

        # Trying new get substruct match method
        product_to_reactant_mapping: Dict[int,int] = self.use_fmcs(reactant)

        attachment_idxs_list = []
        for bond in self.product.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            i_inSub = i in product_to_reactant_mapping.keys()
            j_inSub = j in product_to_reactant_mapping.keys()
            if int(i_inSub) + int(j_inSub) == 1:
                if i_inSub:
                    product_idx = i # This is the attachment point within product
                    try:
                        reactant_idx = product_to_reactant_mapping[i]
                    except (ValueError, KeyError):
                        reactant_idx = None
                else:
                    product_idx = j # This is the attachment point within product
                    try:
                        reactant_idx = product_to_reactant_mapping[j]
                    except (ValueError, KeyError):
                        reactant_idx = None
                attachment_idxs_list.append((product_idx, reactant_idx))

        if len(attachment_idxs_list) == 0 and self.reaction_name == 'Amide_Schotten-Baumann_with_amine':  # run if no other attachments found
            # replace halide with dummy atom
            reactant = self._replace_halide_with_dummy(reactant)
            dummy_idx_list, neig_idx_list = self._find_attachment_id_from_dummy(reactant, dummy_symbol=DUMMY_SYMBOL)
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == 'Amidation':  # run if no other attachments found
            # replace hydroxyl in carboxylic acid with dummy atom
            reactant = self._replace_carboxylic_acid_hydroxy_with_dummy(reactant)
            dummy_idx_list, neig_idx_list = self._find_attachment_id_from_dummy(reactant, dummy_symbol=DUMMY_SYMBOL)
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)':
            # it is probably having trouble on the sulfonyl halide group
            # replace halide with dummy atom
            reactant = self._replace_halide_with_dummy(reactant, atom_to_check_by_symbol='S')
            dummy_idx_list, neig_idx_list = self.find_attachment_id_from_dummy(reactant, dummy_symbol=DUMMY_SYMBOL)
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == "Buchwald-Hartwig_amination":
            # replace halide with dummy atom
            reactant = self._replace_halide_with_dummy(reactant)
            dummy_idx_list, neig_idx_list = self._find_attachment_id_from_dummy(reactant, dummy_symbol=DUMMY_SYMBOL)
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == "Epoxide_+_amine_coupling":
            # Find epoxide group, then find which carbon it is connected to. That is the attachment point.
            # This is a hacky solution, but it works for now.
            epoxide_pattern = Chem.MolFromSmarts("[O]1[C][C]1")
            matches = reactant.GetSubstructMatches(epoxide_pattern)
            if matches:
                for match in matches:
                    atoms_to_remove = list(match)
                    for atom in atoms_to_remove:
                        if reactant.GetAtomWithIdx(atom).GetSymbol() == 'C':
                            if reactant.GetAtomWithIdx(
                                    atom).GetTotalNumHs() == 1:  # Find carbon with one hydrogen, not on end
                                attach_idx = atom
                                return [attach_idx]

        exact_attachment = [x[1] for x in attachment_idxs_list]

        return exact_attachment

    def find_attachment_ids_for_all_reactants(self) -> Dict[Chem.Mol, List[Tuple[int, int]]]:
        """
        This function is used to find the attachment indices of all reactants in the reaction.
        :returns a list of lists, each containing tuples of attachment indices for each reactant
        """
        all_attachments = {}
        for reactant in self.reactants:
            attachments: List[int] = self.find_attachment_id_for_reactant(reactant)
            all_attachments[reactant] = attachments
        self.all_attach_ids = all_attachments
        if any(len(attach_ids) == 0 for attach_ids in self.all_attach_ids.values()):
            raise ReactionError("No attachment points found for reaction {}".format(self.reaction_name))
            return None
        return self.all_attach_ids

    def find_reaction_atoms_for_all_reactants(self) -> List[Dict]:
        """
        This function is used to find the reaction atoms of both reactants. And how those atoms correspond to the SMARTS
        pattern associated with the reaction.
        """
        # check reactant smarts in both reactants
        matched_reactants: Dict[str, Tuple[Chem.Mol, List[int]]] = self.smarts_handler.assign_reactants_w_rxn_smarts(self.all_attach_ids, self.reaction_name)
        self.matched_smarts_to_reactant = matched_reactants
        return self.matched_smarts_to_reactant