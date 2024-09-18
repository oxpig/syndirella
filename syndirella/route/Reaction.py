#!/usr/bin/env python3
"""
syndirella.route.Reaction.py

This module contains the Reaction class. One instance of this object is used to describe a single reaction.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS
from typing import (List, Dict, Tuple)
from syndirella.SMARTSHandler import SMARTSHandler
from syndirella.error import ReactionError
import syndirella.fairy as fairy
import logging


class Reaction():
    """
    This class contains information about the reaction. It is used to find the atoms on the reactants involved in
    the reaction, check the reaction atoms are the ones that are supposed to be connected.
    """

    def __init__(self,
                 product: Chem.Mol,
                 reactants: List[Chem.Mol],
                 reaction_name: str,
                 smarts_handler: SMARTSHandler,
                 route_uuid: str):
        self.route_uuid: str = route_uuid
        self.logger = logging.getLogger(f"{__name__}")
        self.scaffold: Chem.Mol = product
        self.reactants: List[Chem.Mol] = reactants
        self.reaction_name: str = reaction_name
        self.smarts_handler: SMARTSHandler = smarts_handler
        self.reaction_pattern: Chem.rdChemReactions = self.smarts_handler.reaction_smarts[self.reaction_name]
        self.all_attach_ids: Dict[Chem.Mol, List[int]] | None = None
        self.matched_smarts_to_reactant: Dict[str, Tuple[Chem.Mol, Tuple[int], str]] | None = None
        self.matched_smarts_index_to_reactant: Dict[int, Tuple[Chem.Mol, Tuple[int], str]] | None = None
        self.additional_rxn_options: List[Dict[str, str]] = fairy.load_additional_rxn_options()

    def alt_reaction(self) -> List[Dict[str, str]] | None:
        """
        This function is used to determine if an additional reaction is specified.
        """
        rxns = []
        for rxn in self.additional_rxn_options:
            if rxn['name'] == self.reaction_name:
                rxns.append(rxn)
        return rxns

    def get_additional_reactions(self) -> List[Tuple[str, Tuple[str, str]]] | None:
        """
        This function edits the reactants to make an additional reaction.
        """
        alt_rxns: List[Dict[str, str]] = self.alt_reaction()
        new_reactions: List[Tuple[str, Tuple[str, str]]] = []
        for alt_rxn in alt_rxns: # go through all additional reactions
            for r_index in self.matched_smarts_index_to_reactant: # go through each reactant
                if alt_rxn['reactant_id_to_replace'] == r_index:
                    new_reactant: Chem.Mol = self.edit_reactant(self.matched_smarts_index_to_reactant[r_index],
                                                                alt_rxn['reactant_smarts_to_replace_with'],
                                                                alt_rxn['reactant_smarts_to_replace'],
                                                                alt_rxn['replacement_connecting_atom_id'])
                    if new_reactant:
                        other_reactant: Chem.Mol = next(value for key, value in self.matched_smarts_index_to_reactant.items() if key != r_index)[0]
                        # check that new reactants can produce a scaffold
                        reaction_name: str = alt_rxn['replace_with']
                        can_react: bool = self.check_reaction_can_produce_product(new_reactant, other_reactant, reaction_name)
                        if can_react:
                            self.logger.info(f"Additional reaction {reaction_name} for {self.reaction_name} found and "
                                             f"validated.")
                            new_reactant: str = Chem.MolToSmiles(new_reactant)
                            other_reactant: str = Chem.MolToSmiles(other_reactant)
                            reactants: Tuple[str, str] = (new_reactant, other_reactant)
                            new_reactions.append((reaction_name, reactants))
        return new_reactions

    def check_reaction_can_produce_product(self,
                                           new_reactant: Chem.Mol,
                                           other_reactant: Chem.Mol,
                                           reaction_name: str) -> bool:
        """
        This function is used to check if the new reactants can simply produce a scaffold from the labeled reaction.
        """
        reaction: Chem.rdChemReactions = self.smarts_handler.reaction_smarts[reaction_name]
        new_reactant_combo: Tuple[Chem.Mol] = (new_reactant, other_reactant)
        products: Tuple[Chem.Mol] = reaction.RunReactants(new_reactant_combo)
        if len(products[0]) == 0:
            self.logger.error(f"Additional reaction {reaction_name} cannot produce a scaffold from the new reactants "
                              f"{Chem.MolToSmiles(new_reactant)} {Chem.MolToSmiles(other_reactant)}.")
            return False
        return True

    def edit_reactant(self,
                      reactant: Tuple[Chem.Mol, Tuple[int], str],
                      new_reactant_smarts: str,
                      to_replace_smarts: str,
                      replacement_connecting_atom_id: int) -> Chem.Mol | None:
        """
        Directly edits the reactant to replace SMARTS matched atoms to new SMARTS.
        """
        to_replace: Chem.Mol = Chem.MolFromSmarts(to_replace_smarts)
        to_replace_with: Chem.Mol = Chem.MolFromSmarts(new_reactant_smarts)
        reactant_mol: Chem.Mol = reactant[0]
        replaced_reactants: Tuple[Chem.Mol] = Chem.ReplaceSubstructs(mol=reactant_mol,
                                                                     query=to_replace,
                                                                     replacement=to_replace_with,
                                                                     replacementConnectionPoint=replacement_connecting_atom_id,
                                                                     replaceAll=True)
        for replaced_reactant in replaced_reactants:
            # check mol can be sanitized
            if not fairy.check_mol_sanity(replaced_reactant):
                self.logger.error(f"Cannot sanitize new reactant for reaction {self.reaction_name}")
                return None
            # get attach ids
            attach_ids: List[int] = []
            for mol in self.all_attach_ids.keys():
                if fairy.calculate_tanimoto(mol, reactant_mol) == 1.0: # get attach ids of the right reactant replaced
                    attach_ids: List[int] = self.all_attach_ids[mol]
            if len(attach_ids) == 0:
                self.logger.error(f"No attachment points found for reaction {self.reaction_name}")
                return None
            # get atom ids that match new smarts pattern
            matched_atoms: Tuple[int] = replaced_reactant.GetSubstructMatch(to_replace_with)
            # get all atoms that matched atoms are bonded to since SMARTS might not match attachment points
            for atom in matched_atoms:
                for bond in replaced_reactant.GetAtomWithIdx(atom).GetBonds():
                    if bond.GetBeginAtomIdx() in matched_atoms:
                        attach_ids.append(bond.GetEndAtomIdx())
                    elif bond.GetEndAtomIdx() in matched_atoms:
                        attach_ids.append(bond.GetBeginAtomIdx())
            if any(atom in matched_atoms for atom in attach_ids): # If any attach ids are in the SMARTS of the new reactant
                return replaced_reactant
        return None

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
        mcs = rdFMCS.FindMCS([reactant, self.scaffold])
        mcs_smarts = Chem.MolFromSmarts(mcs.smartsString)
        product_matches = self.scaffold.GetSubstructMatches(mcs_smarts)
        # Can return multiple product_matches, only take the first one
        product_match = product_matches[0]
        reactant_matches = reactant.GetSubstructMatches(mcs_smarts)
        reactant_match = reactant_matches[0]
        product_to_reactant_mapping = dict(zip(product_match, reactant_match))
        return product_to_reactant_mapping

    def find_attachment_id_for_reactant(self, reactant: Chem.Mol) -> List[int] | None:
        """
        This function is used to find the attachment indices of a single reactant in the reaction.
        :param reactant: a single reactant molecule
        :returns a list of tuples (attachmentIdx_whole, attachmentIdx_subMol)
        """
        # product_to_reactant_mapping = self.scaffold.GetSubstructMatch(reactant)

        # Trying new get substruct match method
        product_to_reactant_mapping: Dict[int, int] = self.use_fmcs(reactant)

        attachment_idxs_list = []
        for bond in self.scaffold.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            i_inSub = i in product_to_reactant_mapping.keys()
            j_inSub = j in product_to_reactant_mapping.keys()
            if int(i_inSub) + int(j_inSub) == 1:
                if i_inSub:
                    product_idx = i  # This is the attachment point within scaffold
                    try:
                        reactant_idx = product_to_reactant_mapping[i]
                    except (ValueError, KeyError):
                        reactant_idx = None
                else:
                    product_idx = j  # This is the attachment point within scaffold
                    try:
                        reactant_idx = product_to_reactant_mapping[j]
                    except (ValueError, KeyError):
                        reactant_idx = None
                attachment_idxs_list.append((product_idx, reactant_idx))

        if len(attachment_idxs_list) == 0 and self.reaction_name == 'Amide_Schotten-Baumann_with_amine':  # run if no other attachments found
            # replace halide with dummy atom
            reactant = self._replace_halide_with_dummy(reactant)
            dummy_idx_list, neig_idx_list = self._find_attachment_id_from_dummy(reactant, dummy_symbol="*")
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == 'Amidation':  # run if no other attachments found
            # replace hydroxyl in carboxylic acid with dummy atom
            reactant = self._replace_carboxylic_acid_hydroxy_with_dummy(reactant)
            dummy_idx_list, neig_idx_list = self._find_attachment_id_from_dummy(reactant, dummy_symbol="*")
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)':
            # it is probably having trouble on the sulfonyl halide group
            # replace halide with dummy atom
            reactant = self._replace_halide_with_dummy(reactant, atom_to_check_by_symbol='S')
            dummy_idx_list, neig_idx_list = self.find_attachment_id_from_dummy(reactant, dummy_symbol="*")
            print(neig_idx_list)
            print(dummy_idx_list)
            return neig_idx_list[0]

        if len(attachment_idxs_list) == 0 and self.reaction_name == "Buchwald-Hartwig_amination":
            # replace halide with dummy atom
            reactant = self._replace_halide_with_dummy(reactant)
            dummy_idx_list, neig_idx_list = self._find_attachment_id_from_dummy(reactant, dummy_symbol="*")
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

    def find_attachment_ids_for_all_reactants(self):
        """
        This function is used to find the attachment indices of all reactants in the reaction.
        :returns a list of lists, each containing tuples of attachment indices for each reactant
        """
        all_attachments = {}
        for reactant in self.reactants:
            attachments:  List[int] | None = self.find_attachment_id_for_reactant(reactant)
            all_attachments[reactant] = attachments
        if any(len(attach_ids) == 0 for attach_ids in all_attachments.values()):
            raise ReactionError(message=f"No attachment points found for reaction {self.reaction_name}",
                                mol=self.scaffold,
                                route_uuid=self.route_uuid)
        self.all_attach_ids: Dict[Chem.Mol, List[int]] = all_attachments

    def format_matched_smarts_to_index(self, matched_reactants: Dict[str, Tuple[Chem.Mol, List[int], str]]) -> Dict[int,
    Tuple[Chem.Mol, List[int], str]] | None:
        """
        Formats matched smarts to reactant by using reactant index in smarts as the key.
        """
        matched_smarts_index_to_reactant = {}
        for smarts, mol_id_rname in matched_reactants.items():
            mol, ids, _ = mol_id_rname
            if '1' in mol_id_rname[-1]:
                matched_smarts_index_to_reactant[1] = (mol, ids, smarts)
            elif '2' in mol_id_rname[-1]:
                matched_smarts_index_to_reactant[2] = (mol, ids, smarts)
            else:
                self.logger.error(f"Cannot correctly store matched smarts to reactant index for {self.reaction_name}.")
                return None
        return matched_smarts_index_to_reactant


    def find_reaction_atoms_for_all_reactants(self):
        """
        This function is used to find the reaction atoms of both reactants. And how those atoms correspond to the SMARTS
        pattern associated with the reaction.
        """
        # check reactant smarts in both reactants
        matched_reactants: Dict[str, Tuple[Chem.Mol, List[int], str]] | None = (
            self.smarts_handler.assign_reactants_w_rxn_smarts(product=self.scaffold,
                                                              reactant_attach_ids=self.all_attach_ids,
                                                              reaction_name=self.reaction_name))
        self.matched_smarts_index_to_reactant: Dict[int, Tuple[Chem.Mol, List[int], str]] = (
            self.format_matched_smarts_to_index(matched_reactants))
        self.matched_smarts_to_reactant = matched_reactants
        if len(self.matched_smarts_to_reactant) == 0:
            raise ReactionError(message=f"No reaction atoms found for reaction {self.reaction_name}",
                                mol=self.scaffold,
                                route_uuid=self.route_uuid)
