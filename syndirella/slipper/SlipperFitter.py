#!venv/bin/env python3
"""
slipper_fitter/CobblersWorkshop.py

This module contains the SlipperFitter class.
"""
from typing import (List, Dict, Tuple, Union, Optional)
from fragmenstein import Laboratory, Wictor, Igor
import pandas as pd
from fragmenstein.laboratory.validator import place_input_validator
import os, logging
import time
from syndirella.slipper import intra_geometry, flatness
import syndirella.fairy as fairy
import datetime
from syndirella.error import *
import json
import syndirella.slipper._placement_data as placement_data


class SlipperFitter:
    """
    This class is instantiated to place all the products in the template using hits_path
    """

    def __init__(self,
                 template_path: str,
                 hits_path: str,
                 hits_names: List[str],
                 output_dir: str,
                 id: str | None = None,
                 route_uuid: str | None = None,
                 scaffold_placements: Dict[Chem.Mol, str | None] = None):
        self.template_path: str = template_path
        self.hits_path: str = hits_path
        self.hits_names: List[str] = hits_names
        self.output_dir: str = output_dir
        self.id: str | None = id
        self.route_uuid: str | None = route_uuid
        self.scaffold_placements: Dict[Chem.Mol, str | None] = scaffold_placements

        self.atom_diff_min: int = 0
        self.atom_diff_max: int = 10

        self.final_products_csv_path: str | None = None
        self.final_products_pkl_path: str | None = None
        self.batch_num: int | None = None
        self.final_products: pd.DataFrame | None = None
        self.placements: pd.DataFrame | None = None
        self.merged_placements: pd.DataFrame | None = None
        self.output_path = None
        # Placements variables set
        self.n_cores: int = 8
        self.timeout: int = 240
        self.logger = logging.getLogger(f"{__name__}")

        # variables for output
        self.num_placed: int | None = None
        self.num_successful: int | None = None


    def fit_products(self):
        """
        Main entry to the SlipperFitter class. This function is used to fit the products to the final library.
        """
        self.logger.info('Fitting products, warning: it is expected that the scaffold compound has already passed '
                         'minimisation and intramolecular checks.')
        input_df: pd.DataFrame = self.prep_products()
        placements = self.place_products(input_df)
        self.placements = self.check_intra_geometry(placements)
        self.format_placements()
        self.merged_placements = self.merge_placements()
        return self.merged_placements

    def check_scaffold(self,
                       scaffold: Chem.Mol,
                       scaffold_name: str,
                       scaffold_place_num: int) -> str | None:
        """
        Checks if the scaffold can be minimised (no specific stereoisomer) and passes intermolecular checks.
        If it cannot be minimised after 3 attempts, returns False.
        """
        input_df: pd.DataFrame = self._prep_scaffold_input_df(scaffold=scaffold, scaffold_name=scaffold_name)
        id = fairy.generate_inchi_ID(Chem.MolToSmiles(scaffold, isomericSmiles=False))
        output_path: str = os.path.join(self.output_dir, f'{id}-scaffold-check')
        lab: Laboratory = self.setup_Fragmenstein(output_path)
        for attempt in range(1, scaffold_place_num + 1):
            scaffold_placed: Chem.Mol = self._place_scaffold(lab, input_df)  # None if not minimised
            if scaffold_placed is not None:
                paths = [os.path.join(output_path, 'output', scaffold_name, f'{scaffold_name}.minimised.mol'),
                         os.path.join(output_path, scaffold_name, f'{scaffold_name}.minimised.mol')] # could be two output paths...
                path_exists = [os.path.exists(path) for path in paths]
                if not any(path_exists):
                    self.logger.info(f'Scaffold could not be minimised. Attempt {attempt} of {scaffold_place_num}.')
                    continue
                geometries: Dict = intra_geometry.check_geometry(scaffold_placed,
                                                                 threshold_clash=0.4) # increasing threshold for internal clash
                flat_results: Dict = flatness.check_flatness(scaffold_placed)
                if self._check_intra_geom_flatness_results(geometries=geometries, flat_results=flat_results):
                    self.logger.info(f'Scaffold minimised and passed intramolecular checks.')
                    return paths[0] if path_exists[0] else paths[1]
                else:
                    self.logger.info(f'Scaffold could not pass intramolecular checks. Attempt {attempt} of {scaffold_place_num}.')
            else:
                self.logger.info(f'Scaffold could not be minimised. Attempt {attempt} of {scaffold_place_num}.')
        return None

    def _check_intra_geom_flatness_results(self,
                                           geometries: Dict,
                                           flat_results: Dict) -> bool:
        """
        This function checks the intramolecular geometry and flatness results and returns True if all are True.
        """
        if not geometries['results']['bond_lengths_within_bounds']:
            #self.logger.info('Mol could not pass bond length checks.')
            return False
        if not geometries['results']['bond_angles_within_bounds']:
            #self.logger.info('Mol could not pass bond angle checks.')
            return False
        if not geometries['results']['no_internal_clash']:
            #self.logger.info('Mol could not pass internal clash checks.')
            return False
        if not flat_results['results']['flatness_passes']:
            #self.logger.info('Mol could not pass flatness checks.')
            return False
        return True

    def _prep_scaffold_input_df(self,
                                scaffold: Chem.Mol,
                                scaffold_name: str) -> pd.DataFrame:
        """
        Prepares input dataframe to Fragmenstein for the scaffold compound.
        """
        scaffold_df: pd.DataFrame = pd.DataFrame({'name': [scaffold_name],
                                                  'smiles': [Chem.MolToSmiles(scaffold)]})
        scaffold_df = self.add_hits(scaffold_df)
        scaffold_df = scaffold_df.astype({'name': str, 'smiles': str})
        return scaffold_df

    def _get_scaffold(self, input_df: pd.DataFrame) -> Chem.Mol:
        """
        Get scaffold compound as mol object from the input_df.
        """
        # get flat scaffold smiles
        scaffold_smiles: str = input_df.loc['scaffold' in input_df['name'], 'smiles'].values[0]
        scaffold_mol: Chem.Mol = Chem.MolFromSmiles(scaffold_smiles)
        return scaffold_mol

    def prep_products(self) -> pd.DataFrame:
        """
        This function is used to prepare the inputs for Fragmenstein.
        """
        # input_df is a copy of final_products but removing duplicates of names and other metadata from synthesizer step
        input_df: pd.DataFrame = self._prepare_input_df()
        input_df = self.add_hits(input_df)
        # add typing to columns
        input_df = input_df.astype({'name': str, 'smiles': str})
        # Check if there are any duplicates in the input_df
        if input_df.duplicated(subset='name').any():
            self.logger.critical('There are duplicates of names in the scaffold dataframe to place.')
            raise PlacementError(message='There are duplicates of names in the scaffold dataframe to place.',
                                 inchi=self.id,
                                 route_uuid=self.route_uuid)
        return input_df

    def _prepare_input_df(self) -> pd.DataFrame:
        """
        Creates the input dataframe for the placements. Will remove duplicates of names and other metadata from
        synthesizer step.
        """
        final_products: pd.DataFrame = self.final_products.copy(deep=True)
        final_products = final_products.drop_duplicates(subset='name') # drop duplicates of names first
        self.logger.info(f"Total number of products with enumerated stereoisomers (unique name) found before cutting: {len(final_products)}")
        # cut rows with number of atoms difference greater than atom_diff_max
        self.logger.info(f'Cutting products with number of atoms difference greater than {self.atom_diff_max} and '
                         f'below {self.atom_diff_min} to scaffold.')
        input_df = final_products[
            (final_products['num_atom_diff'] <= self.atom_diff_max) &
            (final_products['num_atom_diff'] >= self.atom_diff_min)]
        self._print_diff(orig_df=final_products, input_df=input_df, verb='Kept')
        # drop columns that are not needed
        input_df = input_df[['name', 'smiles']]
        # place number in batch
        if len(input_df) > self.batch_num:
            self.logger.info(f'Even after cutting products, the number of products is over {self.batch_num}.'
                  f' Only placing the top {self.batch_num} products.')
            input_df = input_df.iloc[:self.batch_num]
        self._print_diff(final_products, input_df, verb='Placing')
        return input_df

    def _print_diff(self,
                    orig_df: pd.DataFrame,
                    input_df: pd.DataFrame,
                    verb: str = None):
        """
        This function is used to print the difference between the original number of analogues and the number of
        valid analogues.
        """
        if len(input_df) > len(orig_df):
            self.logger.error("Problem with finding unique analogues. There are more than were in the original list of "
                              "analogues.")
        percent = round(((len(input_df) / len(orig_df)) * 100), 2)
        self.logger.info(f'{verb} {len(input_df)} ({percent}%) enumerated analogue stereoisomers out of {len(orig_df)} '
                         f'enumerated analogue stereoisomers.')

    def add_hits(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """
        This function adds:
            hits_path as mol objects to input_df['hits']
        """
        # load hits_path either from mol or sdf
        if os.path.splitext(self.hits_path)[1] == '.mol':
            self.logger.info('Loading hits from a mol file.')
            hits: List[Chem.Mol] = [Chem.MolFromMolFile(self.hits_path.strip())] # single hit
        else:
            with Chem.SDMolSupplier(self.hits_path.strip()) as sd:
                hits: List[Chem.Mol] = list(sd) # many hits
        # only get hits that exactly match the hit_name in the hits_names
        hits = [hit for hit in hits for hit_name in self.hits_names if hit.GetProp('_Name') == hit_name]
        input_df['hits'] = input_df.apply(lambda row: hits, axis=1)
        return input_df

    def _place_scaffold(self,
                    lab: Laboratory,
                    input_df: pd.DataFrame) -> Chem.Mol:
        """
        Places the scaffold compound, returns the mol object of the scaffold compound if successful else None.
        """
        placements: pd.DataFrame = lab.place(place_input_validator(input_df), n_cores=self.n_cores,
                                             timeout=self.timeout)
        try:
            scaffold: Chem.Mol = placements.at[0, 'minimized_mol']
        except KeyError:
            self.logger.critical('Scaffold could not be minimised.')
            raise ScaffoldPlacementError(message='Scaffold could not be minimised.',
                                         inchi=self.id,
                                         route_uuid=self.route_uuid)
        return scaffold

    def place_products(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """
        This function places products using Fragmenstein.
        """
        start_time = time.time()  # Start timing
        # set up Wictor
        self.output_path: str = os.path.join(self.output_dir, 'output')
        lab: Laboratory = self.setup_Fragmenstein(self.output_path)
        placements: pd.DataFrame = lab.place(place_input_validator(input_df),
                                             n_cores=self.n_cores,
                                             timeout=self.timeout)
        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time  # Calculate elapsed time
        self.logger.info(f"Placing {len(input_df)} run time: {datetime.timedelta(seconds=elapsed_time)}")
        self.num_placed = len(input_df)
        return placements

    def check_intra_geometry(self,
                             placements: pd.DataFrame) -> pd.DataFrame:
        """
        Checks the intramolecular geometry and flatness of double bonds and aromatic rings
        of each mol object using PoseBuster's implemented checks.

        Checks:
        Bond lengths
            The bond lengths in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry
        Bond angles
        	The angles in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry
        Planar aromatic rings
            All atoms in aromatic rings with 5 or 6 members are within 0.25 Å of the closest shared plane
        Planar double bonds
            The two carbons of aliphatic carbon–carbon double bonds and their four neighbours are within 0.25 Å of the closest shared plane
        Internal steric clash
            The interatomic distance between pairs of non-covalently bound atoms is above 0.7 of the lower bound determined by distance geometry
        """
        # get list of mol objects
        mols: List[Chem.Mol] = placements['minimized_mol'].tolist()
        if len(mols) != len(placements):
            self.logger.warning("There are missing mol objects in the placements.")
        intra_geometry_results: List[bool] = []
        flatness_results: List[bool] = []
        for mol in mols:
            if mol is None:
                intra_geometry_results.append(None)
                flatness_results.append(None)
            else:
                geometries: Dict = intra_geometry.check_geometry(mol, threshold_clash=0.4)  # increasing threshold for internal clash
                flat_results: Dict = flatness.check_flatness(mol)
                intra_geometry_results.append(self._check_intra_geom_flatness_results(geometries=geometries, flat_results=flat_results))
        # get list of pass output from intra_geometry for each mol
        placements['intra_geometry_pass'] = intra_geometry_results
        return placements

    def setup_Fragmenstein(self,
                           output_path: str) -> Laboratory:
        """
        This function sets up Fragmenstein to run.
        """
        # Using Wictor to place just by RDKit minimisation
        # make output directory
        os.makedirs(output_path, exist_ok=True)
        Wictor.work_path = output_path
        os.chdir(output_path)  # this does the trick
        Wictor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
        Wictor.monster_joining_cutoff = 5  # Å
        Wictor.quick_reanimation = False  # for the impatient
        Wictor.error_to_catch = Exception  # stop the whole laboratory otherwise
        Wictor.enable_stdout(logging.CRITICAL)
        # Wictor.enable_logfile(os.path.join(self.output_path, f'fragmenstein.log'), logging.ERROR)
        Laboratory.Victor = Wictor
        Igor.init_pyrosetta()
        with open(self.template_path) as fh:
            pdbblock: str = fh.read()
        lab = Laboratory(pdbblock=pdbblock,
                         covalent_resi=None,
                         run_plip=False)
        return lab

    def format_placements(self):
        """
        This function edits the placements by changing outcome column, adding paths to mol files, formatting scaffold
        placements and checking intramolecular geometry.
        """
        self.placements.loc[(self.placements['∆∆G'] < 0) & (self.placements['comRMSD'] < 2), 'outcome'] = 'acceptable'
        self.placements['unminimized_mol'] = self.placements.unminimized_mol.fillna(Chem.Mol())
        self.fix_intxns()
        self.placements = self.add_paths_to_placements(self.placements)
        self.placements = self.format_scaffold_placements(self.placements)
        if self.placements.outcome.value_counts().get('acceptable') is None:
            num_success = 0
            percent_success = 0
        else:
            filtered_df = self.placements[self.placements['intra_geometry_pass'] == True]
            num_success = filtered_df['outcome'].value_counts().get('acceptable', 0)
            percent_success = round((self.placements.outcome.value_counts()['acceptable'] / len(self.placements) * 100)
                                    , 2)
        self.logger.info(f'{num_success} ({percent_success}%) successful placements '
              f'where ∆∆G < 0 and RMSD < 2 Å and passed intramolecular geometry checks.')
        self.num_successful = num_success

    def fix_intxns(self):
        intxn_names: List[tuple] = [c for c in self.placements.columns if isinstance(c, tuple)]
        for intxn_name in intxn_names:
            self.placements[intxn_name] = self.placements[intxn_name].fillna(0).astype(int)

    def format_scaffold_placements(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        This function checks if the scaffold placement was successful. If not, attempt to replace it with the scaffold
        placement performed for the scaffold check. If that placement failed, don't replace.
        """
        placed_scaffolds = df[df['name'].str.contains('scaffold')]
        for i, row in df.iterrows():
            if row['name'] in placed_scaffolds['name'].values:
                # replace with scaffold check placement if scaffold placement failed in this run
                if row['outcome'] != 'acceptable' or row['intra_geometry_pass'] is False or row['minimized_mol'] is None:
                    self.logger.info(f'Placement for {row["name"]} was not successful.')
                    # match which scaffold was placed to final products
                    placed_name: str = row['name']
                    scaffold_mol: Chem.Mol = Chem.MolFromSmiles(row['smiles'])
                    scaffold_path = None
                    for scaffold_check_mol, scaffold_path in self.scaffold_placements.items():
                        if fairy.generate_inchi_ID(mol=scaffold_mol) == fairy.generate_inchi_ID(mol=scaffold_check_mol):
                            scaffold_path = scaffold_path
                            break
                    if scaffold_path is None:
                        self.logger.warning(f'Correct scaffold for {placed_name} could not be found in scaffold-check directory. '
                                            f'Keeping errored scaffold placement from current run.')
                        continue
                    # replace scaffold placement with scaffold check placement
                    ddG, comRMSD, bound, unbound = self.get_scaffold_check_values(scaffold_path)
                    mol = Chem.MolFromMolFile(scaffold_path)
                    geometries: Dict = intra_geometry.check_geometry(mol,
                                                                     threshold_clash=0.4)  # increasing threshold for internal clash
                    flat_results: Dict = flatness.check_flatness(mol)
                    df.at[i, '∆∆G'] = ddG
                    df.at[i, 'comRMSD'] = comRMSD
                    df.at[i, '∆G_bound'] = bound
                    df.at[i, '∆G_unbound'] = unbound
                    if ddG < 0 and comRMSD < 2:
                        df.at[i, 'outcome'] = 'acceptable'
                    else:
                        df.at[i, 'outcome'] = None
                    df.at[i, 'path_to_mol'] = scaffold_path
                    df.at[i, 'intra_geometry_pass'] = self._check_intra_geom_flatness_results(geometries=geometries,
                                                                                              flat_results=flat_results)
                    # remove values from other columns to not mislead stats
                    df.at[i, 'unminimized_mol'] = None
                    df.at[i, 'minimized_mol'] = None
                    df.at[i, 'error'] = None
                    df.at[i, 'runtime'] = None
                    df.at[i, 'LE'] = None
                    self.logger.info(f'Replaced scaffold placement for {placed_name} with scaffold check placement at'
                                     f'{scaffold_path}.')
        return df

    def get_scaffold_check_values(self, scaffold_path: str) -> Tuple[float, float, float, float]:
        """
        This function gets the ∆∆G, comRMSD, ∆G_bound and ∆G_unbound values from the placement done by scaffold check.
        """
        json_file_path = str.replace(scaffold_path, '.mol', '.json')
        if not os.path.exists(json_file_path):
            self.logger.error(f'Could not find json file for scaffold placement in {scaffold_path}.')
            raise PlacementError(message=f'Could not find json file for scaffold placement in {scaffold_path}.',
                                 inchi=self.id,
                                 route_uuid=self.route_uuid)
        with open(json_file_path, 'r') as file:
            data = json.load(file)
        comRMSD = data['mRMSD']
        ddG = placement_data.get_delta_delta_G(data)
        bound, unbound = placement_data.get_bound_unbound(data)
        return ddG, comRMSD, bound, unbound

    def add_paths_to_placements(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds the absolutes paths to the .mol files.
        """
        if self.output_path is None:
            self.logger.error(f'Output path not set for placements of {self.id}.')
            raise PlacementError(message=f'Output path not set for placements of {self.id}.',
                                 inchi=self.id,
                                 route_uuid=self.route_uuid)
        df['paths_to_mols'] = df.apply(lambda x: [os.path.join(self.output_path, x['name'], f'{x["name"]}.minimised.mol'),
                                            os.path.join(self.output_path, 'output', x['name'], f'{x["name"]}.minimised.mol')], axis=1)
        # Add only the first path that exists, or None if neither path exists
        df['path_to_mol'] = df['paths_to_mols'].apply(
            lambda paths: next((path for path in paths if os.path.exists(path)), None))
        df.drop(columns=['paths_to_mols'], inplace=True)
        return df

    def merge_placements(self) -> pd.DataFrame:
        """
        This function merges the fragmenstein output with products csv.
        """
        # Define the columns to drop
        columns_to_drop = ['smiles', 'binary_hits', 'mode', 'disregarded', 'unmin_binary', 'min_binary',
                           'hit_binaries', 'unminimized_mol', 'minimized_mol', 'hit_mols', 'hit_names']
        # Filter out the columns that do not exist in merged_placements
        columns_to_drop = [col for col in columns_to_drop if col in self.placements.columns]
        # Drop the columns if they exist
        self.placements = self.placements.drop(columns=columns_to_drop)
        # merge on name
        merged_placements = pd.merge(self.final_products, self.placements, on='name', how='left')
        return merged_placements

    def save_placements(self, id: str, route_uuid: str):
        """
        This function saves the placements as a .pkl.gz and .csv file.
        """
        self.placements.to_pickle(os.path.join(self.output_dir,
                                               f'{id}_{route_uuid}_fragmenstein_placements.pkl.gz'))
        self.placements.to_csv(os.path.join(self.output_dir,
                                            f'{id}_{route_uuid}_fragmenstein_placements.csv'))
        # get final name same as final products csv
        merged_placements_csv_path = self.final_products_csv_path.split('.')[0] + '_placements.csv'
        merged_placements_pkl_path = self.final_products_pkl_path.split('.')[0] + '_placements.pkl.gz'
        self.merged_placements.to_csv(merged_placements_csv_path)
        self.merged_placements.to_pickle(merged_placements_pkl_path)

