#!venv/bin/env python3
"""
slipper_synthesizer/Labeler.py

This module contains a helper module for the SlipperSynthesizer class to label the products based on labeled atom ids
to expand or not expand.
"""
import os
from rdkit import Chem
from rdkit.Chem import rdFMCS
from typing import (List, Dict, Tuple, Union, Optional)
import pandas as pd
from rdkit.Chem.Draw import rdMolDraw2D
from syndirella.route.Library import Library
import time

class Labeler:
    def __init__(self, products: pd.DataFrame, atom_ids_expansion: Dict[int, bool], library: Library):
        self.products: pd.DataFrame = products
        self.atom_ids_expansion: dict = atom_ids_expansion
        self.library: Library = library
        self.output_dir: str = library.output_dir

    def label_products(self):
        """
        This is the main entry function for the Labeler class.
        """
        self.show_atoms_to_expand_and_not_expand() # save png of colored atoms
        self.show_mcs_on_products() # save png of mcs on products
        self.label_products_with_atom_ids()
        return self.products

    def show_atoms_to_expand_and_not_expand(self):
        """
        This will save a png of the scaffold compound with atoms to expand colored in green and not to expand colored in
        red.
        """
        # get scaffold compound
        base_compound: Chem.Mol = self.library.reaction.scaffold
        # label atom ids
        for i, atom in enumerate(base_compound.GetAtoms()):
            atom_index = atom.GetIdx()
            atom.SetProp("molAtomMapNumber", str(atom_index))
        # assert that atom ids are in the scaffold compound
        assert all([atom_id in [atom.GetIdx() for atom in base_compound.GetAtoms()] for atom_id in
                    self.atom_ids_expansion.keys()]), "Atom ids to expand are not in the scaffold compound."
        atom_ids: List[int] = [atom_id for atom_id in self.atom_ids_expansion.keys() if
                               self.atom_ids_expansion[atom_id] is True or self.atom_ids_expansion[atom_id] is False]
        green = (0, 1, 0)
        red = (1, 0, 0)
        atom_colors: Dict[int, str] = {atom_id: green if self.atom_ids_expansion[atom_id] is True
                                       else red for atom_id in atom_ids}
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(base_compound,
                            highlightAtoms=atom_ids,
                            highlightAtomColors=atom_colors)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        with open(os.path.join(self.output_dir, "base_expansion.svg"), "w") as f:
            f.write(svg)

    def label_products_with_atom_ids(self):
        """
        This function will label the products with the atom ids to expand and not expand.
        """
        # get the atom ids to expand and not expand
        atom_ids_to_expand: List[int] = [atom_id for atom_id in self.atom_ids_expansion.keys() if
                                         self.atom_ids_expansion[atom_id] is True]
        atom_ids_to_not_expand: List[int] = [atom_id for atom_id in self.atom_ids_expansion.keys() if
                                             self.atom_ids_expansion[atom_id] is False]
        # time how long this takes
        start = time.time()
        self.products = self.products.apply(lambda row:
                                            self._label_product_with_atom_ids(atom_ids_to_expand,
                                                                              atom_ids_to_not_expand,
                                                                              row), axis=1)
        end = time.time()
        elapsed_time = end - start
        # Convert elapsed time to hours, minutes, and seconds
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = elapsed_time % 60
        print(f"The function took {hours} hours, {minutes} minutes, and {seconds:.2f} seconds to complete.")

    def _label_product_with_atom_ids(self, atom_ids_to_expand: List[int], atom_ids_to_not_expand: List[int],
                                     row: pd.Series) -> pd.Series:
        """
        This function goes through each scaffold and checks if the compound has been expanded or not.
        """
        results_atom_ids_to_expand: Dict[int, bool] = {}
        results_atom_ids_to_not_expand: Dict[int, bool] = {}
        for atom_id in atom_ids_to_expand:
            results_atom_ids_to_expand[atom_id]: List[bool or str] = self._has_non_mcs_bond(row["smiles"],
                                                                                            self.library.reaction.scaffold,
                                                                                            atom_id)
            row_name = f"expanded_on_atom_{atom_id}"
            row[row_name] = results_atom_ids_to_expand[atom_id]
        for atom_id in atom_ids_to_not_expand:
            results_atom_ids_to_not_expand[atom_id]: List[bool or str] = self._has_non_mcs_bond(row["smiles"],
                                                                                                self.library.reaction.scaffold,
                                                                                                atom_id)
            row_name = f"expanded_on_atom_{atom_id}"
            row[row_name] = results_atom_ids_to_not_expand[atom_id]
        # Qualify if good expansion on if only expanding on good atoms and not expanding on bad atoms
        import math

        # Assuming results_atom_ids_to_not_expand is a dictionary
        # and row is a dictionary representing a row in a dataframe or similar structure

        if any(value == 'non_mcs_match' for value in results_atom_ids_to_expand.values()) or \
                any(value == 'non_mcs_match' for value in results_atom_ids_to_not_expand.values()):
            row["bad_expansion"] = 'non_mcs_match'
            row["good_expansion"] = 'non_mcs_match'
        elif any(value for value in results_atom_ids_to_expand.values()) and \
            not any(value for value in results_atom_ids_to_not_expand.values()):
            # if expanding on any good atom --> good, if not expanding on any bad atom --> good
            row["bad_expansion"] = False
            row["good_expansion"] = True
        elif any(value for value in results_atom_ids_to_not_expand.values()): # if expanding on any bad atom --> bad
            row["bad_expansion"] = True
            row["good_expansion"] = False
        else: # if there is expansion on both bad and good --> indeterminate --> None
            row["bad_expansion"] = None
            row["good_expansion"] = None
        return row

    def _has_non_mcs_bond(self, mol_smiles: str, ref_mol: Chem.Mol, atom_id_ref: int) -> bool or str:
        """
        Check if the specified atom in the given molecule has a bond that is not in the MCS with a reference molecule.

        :param mol_smiles: SMILES string of the molecule to check
        :param ref_smiles: SMILES string of the reference molecule
        :param atom_id: Atom ID to check in the given molecule
        :return: True if the atom has a non-MCS bond, False otherwise
        """
        mol = Chem.MolFromSmiles(mol_smiles)
        # Find MCS between mol and ref_mol
        mcs_res = rdFMCS.FindMCS([mol, ref_mol])
        mcs_mol = Chem.MolFromSmarts(mcs_res.smartsString)
        # check that mcs matches the ref_mol exactly so we are accurately comparing atom indicies
        if mcs_mol.GetNumAtoms() != ref_mol.GetNumAtoms():
            return 'non_mcs_match'
        # Get atom mapping for the MCS in ref
        ref_match = {ref_idx: mcs_idx for mcs_idx, ref_idx in enumerate(ref_mol.GetSubstructMatch(mcs_mol))}
        # Get atom mapping for the MCS in mol
        mcs_match = mol.GetSubstructMatch(mcs_mol)
        # if atom not found in ref_match, ignore
        # Find mcs atom id
        try:
            mcs_atom_id = ref_match[atom_id_ref]
        except KeyError:
            return 'non_mcs_match' # corresponds to not found
        # Find atom id in mol to check from mcs atom id
        atom_id = mcs_match[mcs_atom_id]
        # Check if the atom_id is valid
        if atom_id >= mol.GetNumAtoms():
            raise ValueError("Invalid atom ID.")
        # Check bonds of the specified atom
        for bond in mol.GetAtomWithIdx(atom_id).GetBonds():
            # Get the indices of the bonded atoms
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            # Check if this bond is in the MCS
            if begin_atom_idx not in mcs_match or end_atom_idx not in mcs_match:
                # Bond is not in MCS
                return True
        # No non-MCS bonds found
        return False

    def show_mcs_on_products(self):
        """
        This function will output a png of the MCS on the products.
        """
        return NotImplementedError
