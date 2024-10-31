#!venv/bin/env python3
"""
syndirella.slipper.Slipper.py

This module contains the Slipper class. A slipper in this metaphor is the set of molecules that is the
scaffold of a reaction.
"""
from typing import List, Dict, Optional
import pandas as pd

from syndirella.error import NoScaffold, NoToHippo
from syndirella.slipper.slipper_synthesizer.SlipperSynthesizer import SlipperSynthesizer
from syndirella.route.Library import Library
from syndirella.slipper.SlipperFitter import SlipperFitter
from syndirella.slipper._placement_data import get_placement_data
import os, shutil
import glob
from rdkit import Chem
from rdkit.Chem import inchi
import logging


class Slipper:
    """
    This class is instantiated to represent all products for a step in a route.
    """

    def __init__(self,
                 *,
                 library: Library,
                 template: str = None,
                 hits_path: str = None,
                 hits_names: List[str] = None,
                 batch_num: int = None,
                 atoms_ids_expansion: dict = None,
                 additional_info: dict = None,
                 scaffold_placements: Dict[Chem.Mol, str | None] = None):
        self.products: pd.DataFrame = None
        self.library: Library = library
        self.route_uuid: str = library.route_uuid
        self.output_dir: str = library.output_dir
        self.final_products_pkl_path: str = None
        self.final_products_csv_path: str = None
        self.scaffold_placements: Dict[Chem.Mol, str | None] = scaffold_placements
        # need Fragmenstein information
        self.template: str = template  # path to pdb file
        self.hits_path: str = hits_path  # path to .sdf or .mol file
        self.hits_names: List[str] = hits_names  # name of fragments
        self.batch_num: int = batch_num
        self.atoms_ids_expansion: dict = atoms_ids_expansion
        self.placements: pd.DataFrame = None
        self.output_path: str = None
        self.additional_info: dict = additional_info
        self.logger = logging.getLogger(f"{__name__}")

        # stats for output
        self.num_placed: int | None = None
        self.num_successful: int | None = None
        self.to_hippo_path: str | None = None
        self.num_unique_products: int | None = None # number of unique products not including stereoisomers found at the end of the route.
        self.num_products_enumstereo: int | None = None # number of products after stereoisomer enumeration found at the end of the route.

    def get_products(self) -> pd.DataFrame and str:
        """
        Main entry to the Slipper class. This function is used to get the products the self.library object.
        """
        slipper_synth = SlipperSynthesizer(self.library,
                                           self.output_dir,
                                           self.atoms_ids_expansion,
                                           self.additional_info)
        self.products: pd.DataFrame = slipper_synth.get_products()
        if self.atoms_ids_expansion is not None:
            slipper_synth.label_products()
        slipper_synth.save_products()
        self.final_products_pkl_path: str = slipper_synth.final_products_pkl_path
        self.final_products_csv_path: str = slipper_synth.final_products_csv_path

        # stats for output
        self.num_unique_products = slipper_synth.num_unique_products
        self.num_products_enumstereo = slipper_synth.num_products_enumstereo

        return self.products

    def place_products(self):
        """
        This function is used to place the products with Fragmenstein.
        """
        slipper_fitter = SlipperFitter(template_path=self.template,
                                       hits_path=self.hits_path,
                                       hits_names=self.hits_names,
                                       output_dir=self.output_dir,
                                       route_uuid=self.route_uuid,
                                       id=self.library.id,
                                       scaffold_placements=self.scaffold_placements)
        slipper_fitter.atom_diff_min = self.library.atom_diff_min
        slipper_fitter.atom_diff_max = self.library.atom_diff_max
        slipper_fitter.final_products = self.products # products with enumerated stereoisomers from final library
        slipper_fitter.batch_num = self.batch_num
        slipper_fitter.final_products_pkl_path = self.final_products_pkl_path
        slipper_fitter.final_products_csv_path = self.final_products_csv_path

        self.placements: pd.DataFrame = slipper_fitter.fit_products()

        slipper_fitter.save_placements(id=self.library.id, route_uuid=self.route_uuid)
        self.output_path: str = slipper_fitter.output_path

        # stats for output
        self.num_placed = slipper_fitter.num_placed
        self.num_successful = slipper_fitter.num_successful

        return self.placements

    def write_products_to_hippo(self) -> str:
        """
        Writes a dataframe that contains the values needed for HIPPO db input.

        Returns:
            path: str : the path to the saved dataframe
        """
        if self.placements is None:
            self.logger.critical("Placements need to be run first before writing HIPPO output.")
            return None
        # cut placements to those that were placed by batch_num
        placements: pd.DataFrame = self.placements.iloc[:self.batch_num]

        hippo_path: str = os.path.join(self.output_dir, f'{self.library.id}_{self.route_uuid}_to_hippo.pkl.gz')
        # get all products dfs in /extra
        products_files: List[str] = glob.glob(f"{self.output_dir}/extra/*{self.route_uuid}*products*.pkl*")
        product_dfs: Dict[int, pd.DataFrame] = self._load_products_dfs(products_files)
        # make HIPPO output dataframe of these specific products
        hippo_df = self._structure_products_for_hippo(placements_df=placements,
                                                      product_dfs=product_dfs)
        # load file if it already exists
        if os.path.isfile(hippo_path):
            previous_hippo = pd.read_pickle(hippo_path)
            hippo_df = pd.concat([previous_hippo, hippo_df], axis=1, ignore_index=True)
        # save the dataframe
        hippo_df.to_pickle(hippo_path)
        self.logger.info(f"Saved HIPPO output to {hippo_path}")
        self.to_hippo_path = hippo_path
        self.check_scaffold_in_hippo(hippo_df, hippo_path)
        return hippo_path

    def check_scaffold_in_hippo(self, hippo_df: pd.DataFrame, hippo_path: str):
        """
        Checks if there is a scaffold in the scaffold names of the HIPPO output.
        """
        if not any('scaffold' in name for name in hippo_df[f'{self.library.num_steps}_product_name']):
            self.logger.warning("Scaffold was not found in the scaffold names of the HIPPO output.")
            raise NoScaffold(message=f"Scaffold was not found in the scaffold names of the HIPPO output at {hippo_path}."
                                     f"Most likely due to an incorrectly written SMIRKS.",
                             route_uuid=self.route_uuid,
                             inchi=self.library.id)

    def _load_products_dfs(self, products_files: List[str]) -> Dict[int, pd.DataFrame]:
        """
        Load the products dataframes from the files in the /extra directory and putting into dict where key is step.
        """
        product_dfs: Dict[int, pd.DataFrame] = {}
        for file in products_files:
            df = pd.read_pickle(file)
            if len(df) == 0:
                self.logger.info(f"Empty dataframe found in {file}. Continuing with next file to structure hippo outputs")
                continue
            step = df['step'].iloc[0]
            product_dfs[step] = df
        if len(product_dfs) != self.library.num_steps - 1:
            found_steps = list(product_dfs.keys())
            missing_steps = [step for step in range(1, self.library.num_steps) if step not in found_steps]
            self.logger.critical("Not all steps have findable products dataframes or non empty dataframes. "
                                 f"Missing steps: {missing_steps}")
            raise NoToHippo(message=f"Not all steps have findable products dataframes or non empty dataframes. "
                                    f"Missing steps: {missing_steps}",
                            route_uuid=self.route_uuid,
                            inchi=self.library.id)
        return product_dfs

    def _structure_products_for_hippo(self,
                                      placements_df: pd.DataFrame,
                                      product_dfs: Dict[int, pd.DataFrame]) -> pd.DataFrame:
        """
        Structures the placements or products df for HIPPO output.
        """
        hippo_dfs: Dict[int, pd.DataFrame] = {}
        for step, products_df in product_dfs.items():
            hippo_df_step: pd.DataFrame = self._structure_step_for_hippo(step, products_df)
            hippo_dfs[step] = hippo_df_step
        hippo_dfs[self.library.num_steps] = self._structure_step_for_hippo(self.library.num_steps, placements_df)
        hippo_df = self._put_hippo_dfs_together(hippo_dfs)
        return hippo_df

    def calculate_inchi_similarity(self,
                                   smiles1: str,
                                   smiles2: str) -> int:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        # Check if either molecule failed to be created from the SMILES
        if mol1 is None or mol2 is None:
            return 0  # Indicating an error or non-comparability
        inchi1 = inchi.MolToInchi(mol1)
        inchi2 = inchi.MolToInchi(mol2)
        # Directly compare the InChI strings for equality
        similarity = 1 if inchi1 == inchi2 else 0
        return similarity

    # Function to find matches
    def find_matches(self,
                     row,  # row of current step
                     step,
                     df_previous_step: pd.DataFrame) -> List[str]:
        product_matches = None  # Store matching scaffold names
        for _, product_row in df_previous_step.iterrows():
            similarity = self.calculate_inchi_similarity(
                row[f'{step + 1}_r{row[f"{step + 1}_r_previous_product"]}_smiles'],
                product_row[f'{step}_product_smiles'])
            if similarity == 1:
                product_matches = product_row[f'{step}_product_name']
                break
        return product_matches

    def _put_hippo_dfs_together(self,
                                hippo_dfs: Dict[int, pd.DataFrame]) -> pd.DataFrame:
        """
        Puts the HIPPO dataframes together by matching on each reaction scaffold to the correct previous step's reactant.
        """
        # get the last step's dataframe
        hippo_df_step_last = hippo_dfs[self.library.num_steps]
        for step in range(1, self.library.num_steps)[::-1]:  # iterate through the steps in reverse
            # get the step's dataframe
            hippo_df_stepx = hippo_dfs[step]
            # Find the matching scaffold names for the reactant in this new step
            hippo_df_step_last[f'{step}_product_name'] = hippo_df_step_last.apply(self.find_matches,
                                                                                  step=step,
                                                                                  df_previous_step=hippo_df_stepx,
                                                                                  axis=1)
            # What happens if there are null scaffold names?... Still keep row with null scaffold name
            # Join on the scaffold names, has to be right merge because we only care about the products from the last step
            result_df = pd.merge(hippo_df_stepx,
                                 hippo_df_step_last,
                                 left_on=f'{step}_product_name',
                                 right_on=f'{step}_product_name',
                                 how='right')
            # Update the last step dataframe
            hippo_df_step_last = result_df
        # add scaffold compound smiles as first column, get from Inchi Key
        base_compound_smiles: str = Chem.MolToSmiles(self.library.reaction.scaffold)
        hippo_df_step_last.insert(0, 'scaffold_smiles', base_compound_smiles)
        return hippo_df_step_last

    def _structure_step_for_hippo(self,
                                  step: int,
                                  products_df: pd.DataFrame) -> pd.DataFrame:
        """
        Structures the products df for HIPPO output.
        """
        # cut down the dataframe to those with num_atom_diff <= 15
        reaction: str = products_df['reaction'].iloc[0]
        r1_smiles: List[str] = products_df['r1_smiles'].tolist()
        r1_is_previous_product: bool = products_df['r1_is_previous_product'].iloc[0]
        try:
            r2_smiles: List[str] = products_df['r2_smiles'].tolist()
            r2_is_previous_product: bool = products_df['r2_is_previous_product'].iloc[0]
        except KeyError:
            r2_smiles: List[str] = [''] * len(products_df)
            r2_is_previous_product: bool = None
        # find which reactants were products of the previous step
        r_previous_product: int | None = self.which_reactant_was_previous_product(r1_is_previous_product,
                                                                           r2_is_previous_product)
        product_smiles: List[str] = products_df['smiles'].tolist()
        product_names: List[str] = products_df['name'].tolist()
        try:
            flags: List[str] = products_df['flag'].tolist()
            # replace nan with None
            flags = [None if pd.isna(flag) else flag for flag in flags]
        except KeyError:
            flags: List[str | None] = [None] * len(products_df)
        # variables for the last step
        if step == self.library.num_steps:
            num_atom_diff: List[int] = products_df['num_atom_diff'].tolist()
            stereoisomer: List[str] = products_df['stereoisomer'].tolist()
            error: List[str] = products_df['error'].tolist()
            delta_delta_G: List[float] = products_df['∆∆G'].tolist()
            delta_G_bound: List[float] = products_df['∆G_bound'].tolist()
            delta_G_unbound: List[float] = products_df['∆G_unbound'].tolist()
            comRMSD: List[float] = products_df['comRMSD'].tolist()
            regarded: List[bool] = products_df['regarded'].tolist()
            intra_geometry_pass: List[bool] = products_df['intra_geometry_pass'].tolist()
            path_to_mol: List[Optional[str]] = products_df['path_to_mol'].tolist()
        # make HIPPO output dataframe
        hippo_df_step = pd.DataFrame({f'{step}_reaction': reaction,
                                      f'{step}_r1_smiles': r1_smiles,
                                      f'{step}_r2_smiles': r2_smiles,
                                      f'{step}_r_previous_product': r_previous_product,
                                      f'{step}_product_smiles': product_smiles,
                                      f'{step}_product_name': product_names,
                                      f'{step}_flag': flags})
        if step == self.library.num_steps:
            hippo_df_step[f'{step}_num_atom_diff'] = num_atom_diff
            hippo_df_step[f'{step}_stereoisomer'] = stereoisomer
            hippo_df_step[f'error'] = error
            hippo_df_step[f'∆∆G'] = delta_delta_G
            hippo_df_step[f'∆G_bound'] = delta_G_bound
            hippo_df_step[f'∆G_unbound'] = delta_G_unbound
            hippo_df_step[f'comRMSD'] = comRMSD
            hippo_df_step[f'regarded'] = regarded
            hippo_df_step[f'path_to_mol'] = path_to_mol
            hippo_df_step['template'] = self.template
            hippo_df_step[f'intra_geometry_pass'] = intra_geometry_pass
        return hippo_df_step

    def which_reactant_was_previous_product(self,
                                            r1_is_previous_product: bool,
                                            r2_is_previous_product: bool) -> int | None:
        """
        Determine which reactant was the scaffold of the previous step.
        """
        if r1_is_previous_product:
            if r2_is_previous_product:
                self.logger.critical("Both reactants cannot be products of the previous step.")
                return None
            return 1
        elif r2_is_previous_product:
            return 2
        # what if both ?? don't throw error still work
        else:
            return None

    def _delete_file_or_directory(self, path):
        """
        Delete a file or directory at the given path.
        """
        try:
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
                #print(f"Deleted file: {path}")
            elif os.path.isdir(path):
                shutil.rmtree(path)
                #print(f"Deleted directory: {path}")
        except Exception as e:
            self.logger.warning('Failed to delete %s. Reason: %s' % (path, e))

    def _should_delete_file(self, file, suffixes_to_keep):
        """
        Determine if a file should be deleted based on its suffix.
        """
        return not any(file.endswith(suffix) for suffix in suffixes_to_keep)

    def clean_up_placements(self):
        """
        This function is used to remove extra files that are generated by Fragmenstein.
        """
        self.logger.info("Cleaning up placement directory...")
        suffixes_to_keep = ['.minimised.json', '.minimised.mol', '.csv', '.pkl', '.pkl.gz']
        for root, dirs, files in os.walk(self.output_path):
            for file in files:
                if self._should_delete_file(file, suffixes_to_keep):
                    file_path = os.path.join(root, file)
                    self._delete_file_or_directory(file_path)

    def get_placements_df(self) -> pd.DataFrame | None:
        """
        This function is used to get the placements dataframe which is the merged df with the products.
        """
        if self.placements is not None:
            return self.placements
        # otherwise you're gonna have to build it from scratch by going through each file in the output dir
        # get products df
        if self.products is None:
            self.logger.critical("You need to set the products df first.")
            return None
        self.placements = get_placement_data(self.products, self.output_path, self.output_dir)
        return self.placements
