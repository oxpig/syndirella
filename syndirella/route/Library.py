#!venv/bin/env python3
"""
syndirella.route.Library.py

This module contains the Library class. This class contains information about the analogue library. It will create the
analogue library from the Reaction object. It will also store the analogue library as a .pkl file.
"""

import glob
import logging
import os
import pickle
from typing import (List, Dict, Tuple, Optional)

import pandas as pd
from pandas import DataFrame
from rdkit import DataStructs
from rdkit.Chem.FilterCatalog import *

import syndirella.utils.fairy as fairy
from syndirella.database.Postera import Postera
from syndirella.database.Arthor import Arthor
from syndirella.database.Hippo import Hippo
from syndirella.utils.error import SMARTSError, NoReactants, APIQueryError
from .Reaction import Reaction
from . import LibraryConfig
from syndirella.constants import DatabaseSearchTool, DEFAULT_DATABASE_SEARCH_TOOL


class Library:
    """
    This class contains information about the analogue library. It will create the analogue library from the Reaction
    object. It will also store the analogue library as a .pkl file.
    """

    def __init__(self,
                 reaction: Reaction,
                 output_dir: str,
                 id: str,
                 num_steps: int,
                 current_step: int,
                 filter: bool,
                 route_uuid: str,
                 atom_diff_min: int,
                 atom_diff_max: int,
                 db_search_tool: DatabaseSearchTool,
                 elab_single_reactant: bool,
                 elab_single_reactant_int: Optional[int] = None,
                 reference_db: str = None):
        # Create configuration object
        self.config = LibraryConfig(
            output_dir=output_dir,
            id=id,
            atom_diff_min=atom_diff_min,
            atom_diff_max=atom_diff_max,
            filter=filter,
            db_search_tool=db_search_tool,
            num_steps=num_steps,
            current_step=current_step,
            route_uuid=route_uuid,
            elab_single_reactant=elab_single_reactant,
            elab_single_reactant_int=elab_single_reactant_int,
            reference_db=reference_db
        )
        
        # Core attributes
        self.reaction: Reaction = reaction
        self.output_dir: str = os.path.join(self.config.output_dir, self.config.id)
        self.extra_dir_path: str = os.path.join(self.output_dir, "extra")
        self.analogues_dataframes: Dict[str: Tuple[pd.DataFrame, Tuple[str, str]]] = {}
        self.database_search: str = 'postera'
        self.reference_db: str = reference_db

        self.logger = logging.getLogger(f"{__name__}")
        self.r1 = None
        self.r2 = None

    @property
    def id(self) -> str:
        return self.config.id
    
    @property
    def num_steps(self) -> int:
        return self.config.num_steps
    
    @property
    def current_step(self) -> int:
        return self.config.current_step
    
    @property
    def filter(self) -> bool:
        return self.config.filter
    
    @property
    def route_uuid(self) -> str:
        return self.config.route_uuid
    
    @property
    def atom_diff_min(self) -> int:
        return self.config.atom_diff_min
    
    @property
    def atom_diff_max(self) -> int:
        return self.config.atom_diff_max
    
    @property
    def elab_single_reactant_int(self) -> int | None:
        return self.config.elab_single_reactant_int
    
    @elab_single_reactant_int.setter
    def elab_single_reactant_int(self, value: int | None):
        self.config.elab_single_reactant_int = value
    
    @property
    def db_search_tool(self) -> DatabaseSearchTool:
        return self.config.db_search_tool

    def create_library(self):
        """
        This function is used to create the analogue library from the Reaction object.
        It handles the case where only one reactant is elaborated based on self.elab_single_reactant_int.
        """
        for index, (key, value) in enumerate(
                self.reaction.matched_smarts_to_reactant.items()):  # can work for 1 and 2 reactants
            reactant: Chem.Mol = value[0]
            reactant_smarts: str = key
            analogue_prefix: str = value[2]

            if self.elab_single_reactant_int is not None and index != self.elab_single_reactant_int:
                # Run only if elab_single_reactant_int is set at matches the index in matched smarts
                self.logger.info(
                    f"Skipping elaboration for reactant {analogue_prefix}. Storing original reactant only.")
                df, analogue_columns = (
                    self.process_analogues([Chem.MolToSmiles(reactant)],
                                           reactant_smarts,
                                           analogue_prefix,
                                           previous_product=False,
                                           original_reactant=reactant))
                self.save_library(df, analogue_prefix)
            else:
                # Process the reactant as usual
                df, analogue_columns = self.process_reactant(reactant, reactant_smarts, analogue_prefix)
            df: pd.DataFrame
            analogue_columns: Tuple[str, str]
            self.analogues_dataframes[analogue_prefix] = (df, analogue_columns)

    def process_reactant(self, reactant: Chem.Mol, reactant_smarts: str, analogue_prefix: str) -> tuple[
        DataFrame, tuple[str, str]]:
        # search for analogues csv if already created. Perform for all reactants.
        reactant_analogues_path, internal_step, previous_product = self.check_analogues_pkl_exists(analogue_prefix,
                                                                                                   reactant)
        analogues: List[str] | None = None
        
        if reactant_analogues_path is not None:
            analogues: List[str] = self.load_library(reactant_analogues_path,
                                                     analogue_prefix,
                                                     internal_step)
        else:
            # Perform database search based on the selected tool
            if self.db_search_tool == DatabaseSearchTool.ARTHOR:
                self.logger.info(f"Using Arthor for database search for {analogue_prefix}")
                database_search = Arthor()
            elif self.db_search_tool == DatabaseSearchTool.MANIFOLD:
                self.logger.info(f"Using Postera/Manifold for database search for {analogue_prefix}")
                database_search = Postera()
            elif self.db_search_tool == DatabaseSearchTool.HIPPO:
                assert self.reference_db
                self.logger.info(f"Using HIPPO for database search for {analogue_prefix}")
                self.logger.info(f"reference_db={self.reference_db}")
                database_search = Hippo(self.reference_db)
            else:
                self.logger.error(f"Database search tool {self.db_search_tool} not found.")
                raise ValueError(f"Database search tool {self.db_search_tool} not found.")
            
            analogues: List[str] | None = database_search.perform_database_search(
                reactant=reactant,
                reaction_name=self.reaction.reaction_name,
                search_type="superstructure"
            )
            
            if analogues is None:
                self.logger.warning(f"Database search failed for {analogue_prefix}.")
                raise APIQueryError(f"Database search of {self.db_search_tool} failed for {analogue_prefix}.")
        
        df, analogue_columns = self.process_analogues(analogues, reactant_smarts, analogue_prefix, previous_product,
                                                        original_reactant=reactant)
        self.save_library(df, analogue_prefix)
        return df, analogue_columns

    def process_analogues(self,
                          analogues: List[str],
                          reactant_smarts: str,
                          analogue_prefix: str,
                          previous_product: bool,
                          original_reactant: Chem.Mol = None) -> pd.DataFrame:
        """
        This function puts list of analogues in dataframe and does SMART checking to check if the analogues contains
        the SMARTS pattern of the original reactant and against all other reactants SMARTS.
        Also detects protecting groups on the original reactant and filters analogues to only keep those
        that contain the same protecting groups.
        """
        analogues_mols: List[Chem.Mol] = [Chem.MolFromSmiles(analogue) for analogue in analogues]
        analogues_mols: List[Chem.Mol] = fairy.remove_chirality(analogues_mols)
        analogues_mols: List[Chem.Mol] = fairy.remove_repeat_mols(analogues_mols)
        self.print_diff(analogues, analogues_mols, analogue_prefix)
        if self.filter:
            analogues: List[str] = self.filter_analogues(analogues, analogue_prefix)

        # Detect protecting groups on the original reactant and filter analogues
        pg_smarts = fairy.load_substructure_to_include()
        required_pg_names: List[str] = []
        if original_reactant is not None:
            required_pg_names = fairy.detect_protecting_groups(original_reactant, pg_smarts)
            if required_pg_names:
                self.logger.info(
                    f"Protecting group(s) detected on original reactant for {analogue_prefix}: "
                    f"{', '.join(required_pg_names)}. Filtering analogues to match.")
                keep, pg_labels = fairy.filter_analogues_on_protecting_groups(
                    analogues_mols, required_pg_names, pg_smarts)
                num_before = len(analogues_mols)
                analogues_mols = [mol for mol, k in zip(analogues_mols, keep) if k]
                pg_labels = [label for label, k in zip(pg_labels, keep) if k]
                self.logger.info(
                    f"Protecting group filter removed {num_before - len(analogues_mols)} of {num_before} "
                    f"analogues for {analogue_prefix}.")
            else:
                self.logger.info(f"No protecting groups detected on original reactant for {analogue_prefix}.")

        # If no PG filtering was done, generate empty labels
        if not required_pg_names:
            _, pg_labels = fairy.filter_analogues_on_protecting_groups(
                analogues_mols, [], pg_smarts)

        reactant_smarts_mol: Chem.Mol = Chem.MolFromSmarts(reactant_smarts)
        contains_smarts_pattern, num_matches = self.check_analogue_contains_smarts_pattern(analogues_mols,
                                                                                           reactant_smarts_mol)
        contains_smarts_pattern: List[bool]
        num_matches: List[int]
        if len(self.reaction.matched_smarts_to_reactant) == 1:
            contains_other_reactant_smarts_pattern = [False for _ in analogues_mols]
            other_reactant_prefix = None
        else:
            contains_other_reactant_smarts_pattern, other_reactant_prefix = (
                self.check_analogue_contains_other_reactant_smarts_pattern(analogues_mols, reactant_smarts))
        contains_other_reactant_smarts_pattern: List[bool]
        other_reactant_prefix: str
        analogues = [Chem.MolToSmiles(analogue) for analogue in analogues_mols]
        analogues_df = (
            pd.DataFrame({f"{analogue_prefix}_smiles": analogues,
                          f"{analogue_prefix}_mol": analogues_mols,
                          f"{analogue_prefix}_{self.reaction.reaction_name}": contains_smarts_pattern,
                          f"{analogue_prefix}_{self.reaction.reaction_name}_num_matches": num_matches,
                          f"{other_reactant_prefix}_{self.reaction.reaction_name}": contains_other_reactant_smarts_pattern,
                          f"{analogue_prefix}_is_previous_product": previous_product,
                          f"{analogue_prefix}_protecting_group": pg_labels}))
        if self.filter:
            analogues_df['is_PAINS_A'] = False
        return analogues_df, (f"{analogue_prefix}_{self.reaction.reaction_name}",
                              f"{other_reactant_prefix}_{self.reaction.reaction_name}")

    def check_analogue_contains_smarts_pattern(self, analogues_mols: List[Chem.Mol], reactant_smarts_mol: Chem.Mol):
        """
        This function is used to check if the analogues contains the original reactant. Will return
        dictionary of boolean values and the number of matches.
        """
        self.logger.info('Checking if analogues contain SMARTS pattern of original reactant...')
        matching = [bool(analogue.GetSubstructMatches(reactant_smarts_mol)) for analogue in analogues_mols]
        num = [len(analogue.GetSubstructMatches(reactant_smarts_mol)) for analogue in analogues_mols]
        if len(matching) != len(analogues_mols):
            self.logger.error("Problem with finding matches.")
            raise SMARTSError(message="Problem with finding SMARTS matches to analogues.",
                              mol=self.reaction.scaffold,
                              route_uuid=self.route_uuid)
        if len(num) != len(analogues_mols):
            self.logger.error("Problem with finding number of matches.")
            raise SMARTSError(message="Problem with finding number of matches.",
                              mol=self.reaction.scaffold,
                              route_uuid=self.route_uuid)
        return matching, num

    def filter_analogues(self, analogues: List[str], analogue_prefix: str) -> List[str]:
        """
        This function is used to filter analogues.
        """
        if not self.filter:
            return analogues
        
        mols = [Chem.MolFromSmiles(analogue) for analogue in analogues]
        mols = [mol for mol in mols if mol is not None]
        
        if len(mols) == 0:
            return analogues
        
        self.print_diff(mols, mols, analogue_prefix)
        return [Chem.MolToSmiles(mol) for mol in mols]

    def print_diff(self, mols: List[Chem.Mol], valid_mols: List[Chem.Mol], analogue_prefix: str):
        """
        This function is used to print the difference between mols and valid_mols.
        """
        if len(mols) != len(valid_mols):
            self.logger.info(f"Filtered {len(mols) - len(valid_mols)} analogues for {analogue_prefix}.")

    def filter_on_substructure_filters(self, mols: List[Chem.Mol], ) -> List[str]:
        """
        This function is used to filter on substructure filters.
        """
        # Implementation for substructure filtering
        return [Chem.MolToSmiles(mol) for mol in mols]

    def check_analogue_contains_other_reactant_smarts_pattern(self, analogues_mols: List[Chem.Mol],
                                                             reactant_smarts: str) -> List[bool]:
        """
        This function is used to check if analogues contain other reactant SMARTS patterns.
        """
        other_reactant_smarts_pattern = \
            [smarts for smarts in self.reaction.matched_smarts_to_reactant.keys() if not smarts == reactant_smarts][0]
        reactant_smarts_mol = Chem.MolFromSmarts(other_reactant_smarts_pattern)
        if reactant_smarts_mol is None:
            self.logger.error(f"Other reactant SMARTS pattern not found for {self.reaction.reaction_name}.")
            raise SMARTSError(message=f"Other reactant SMARTS pattern not found for {self.reaction.reaction_name}.",
                              mol=self.reaction.scaffold,
                              route_uuid=self.route_uuid)
        
        other_reactant_prefix = self.reaction.matched_smarts_to_reactant[other_reactant_smarts_pattern][2]
        other_reactant_smarts_mol: Chem.Mol = Chem.MolFromSmarts(other_reactant_smarts_pattern)
        return [bool(analogue.GetSubstructMatches(other_reactant_smarts_mol)) for analogue in
                analogues_mols], other_reactant_prefix

    def check_analogue_contains_all_smarts_patterns(self, analogues_mols: List[Chem.Mol]) -> List[
        bool] | NotImplementedError:
        """
        This function is used to check if analogues contain all SMARTS patterns.
        """
        raise NotImplementedError("This function is not implemented yet.")

    def check_analogues_pkl_exists(self, analogue_prefix: str, reactant: Chem.Mol) -> str and bool:
        """
        This function is used to check if the analogue library has already been created and saved as a .pkl file.
        """
        internal_step = False
        previous_product = False
        
        if self.current_step != 1 and (self.current_step == self.num_steps or self.current_step < self.num_steps):
            self.logger.info(
                'Since this is an internal or final step looking for the products .pkl from previous step...')
            try:
                pkl_path: str = self._find_products_pkl_name(reactant)
                if pkl_path is not None:
                    internal_step = True
                    previous_product = True
                    return pkl_path, internal_step, previous_product  # returns path to products .pkl.gz since is reactant
            except NoReactants as e:
                # If NoReactants is raised, it means either no product files found or correct product not contained
                # Don't raise immediately, fall back to analogue search
                self.logger.info(f'No suitable products found from previous step for {analogue_prefix}. Falling back to analogue search.')
    
        if self.current_step != 1 and self.current_step < self.num_steps:
            self.logger.info('Looking for analogue library .pkl.gz if already created...')
            pkl_path: str = self._find_analogues_pkl_name(analogue_prefix)
            if pkl_path is not None:
                internal_step = True
                previous_product = False
                return pkl_path, internal_step, previous_product  # returns path to analogue .pkl.gz and is internal step
    
        if self.current_step == 1:
            self.logger.info('Looking for analogue library .pkl.gz if already created...')
            pkl_path: str = self._find_analogues_pkl_name(analogue_prefix)
            if pkl_path is not None:
                internal_step = False
                previous_product = False
                return pkl_path, internal_step, previous_product  # returns path to analogue .pkl.gz and is first step
    
        previous_product = False
        os.makedirs(self.extra_dir_path, exist_ok=True)
        return None, internal_step, previous_product

    def _find_analogues_pkl_name(self, analogue_prefix: str) -> str:
        """
        Checks if the analogue library was already created and saved as a .pkl.gz file. Returns the path to the
        .pkl.gz file.
        """
        search = (f"{self.extra_dir_path}/"
                  f"{self.id}_{self.route_uuid}_{self.reaction.reaction_name}_{analogue_prefix}_"
                  f"{self.current_step}of{self.num_steps}.pkl.gz")
        self.logger.info(f"Looking for analogue pickle at {search}")
        pkl: List[str] = glob.glob(search)
        if len(pkl) == 1:
            self.logger.info(f"Found {pkl[0]} as the analogue library .pkl from previous step.")
            return pkl[0]

    def _find_products_pkl_name(self, reactant: Chem.Mol) -> str:
        """
        This function is used to find the name of the products .pkl file by comparing the reactant to the scaffold with similarity.
        Raises NoReactants if no suitable product is found.
        """
        product_pkls: List[str] = glob.glob(f"{self.extra_dir_path}/"
                                            f"{self.id}_{self.route_uuid}_*products_{self.current_step - 1}of"
                                            f"{self.num_steps}.pkl.gz")
        
        # If no product files found at all, raise NoReactants (steps not connected)
        if not product_pkls:
            self.logger.info(f"Could not find any products .pkl from previous step.")
            raise NoReactants(
                message="Could not find any products .pkl from previous step.",
                route_uuid=self.route_uuid,
                inchi=self.id,
                mol=self.reaction.scaffold
            )
        
        found_file = None
        for i, path in enumerate(product_pkls):
            self.logger.info(f"Found {path} as the potential products .pkl from previous step.")
            found_file = path
            # Find if the reactant is the same as the scaffold
            df = pd.read_pickle(path)
            # Iterate through the top 100 products
            for product in df["smiles"][:100]:
                product_mol = Chem.MolFromSmiles(product)
                # Calculate similarity score
                similarity_score = DataStructs.FingerprintSimilarity(
                    fairy.get_morgan_fingerprint(reactant),
                    fairy.get_morgan_fingerprint(product_mol))
                # If a perfect match is found, return the path
                if similarity_score == 1:
                    self.logger.info(f"Found {path} as the products .pkl from previous step and correct product is contained.")
                    return path
        
        # If we get here, we found product files but none contained the correct product
        # Don't raise NoReactants here, let the calling function handle it
        self.logger.info(f"Found {len(product_pkls)} product files from previous step, but correct product was not contained in any of them.")
        raise NoReactants(
            message=f"Found product files from previous step, but correct product was not contained.",
            route_uuid=self.route_uuid,
            inchi=self.id,
            mol=self.reaction.scaffold
        )

    def save_library(self, df: pd.DataFrame, analogue_prefix: str):
        """
        This function is used to save the analogue library as a .pkl file in self.extra_dir_path
        """
        os.makedirs(f"{self.extra_dir_path}", exist_ok=True)
        pkl_name = f"{self.id}_{self.route_uuid}_{self.reaction.reaction_name}_{analogue_prefix}_{self.current_step}of{self.num_steps}.pkl.gz"
        self.logger.info(f"Saving {analogue_prefix} analogue library to {self.extra_dir_path}/{pkl_name} \n")
        df.to_pickle(f"{self.extra_dir_path}/{pkl_name}")

    def load_library(self,
                     reactant_analogues_path: str,
                     analogue_prefix: str,
                     internal_step: bool) -> List[str]:
        """
        This function is used to load the analogue library from a .pkl file.
        """
        try:
            # find column with analogue smiles
            df = pd.read_pickle(reactant_analogues_path)
            if len(df) == 0:
                self.logger.critical(f"Empty analogue library at {reactant_analogues_path}. Stopping...")
                raise NoReactants(message=f"Empty analogue library at {reactant_analogues_path}",
                                  route_uuid=self.route_uuid,
                                  inchi=self.id,
                                  mol=self.reaction.scaffold)
            if not internal_step:
                analogues = df[f"{analogue_prefix}_smiles"].tolist()
                return analogues
            if internal_step:
                analogues = df["smiles"].tolist()
                return analogues
        except KeyError:
            self.logger.critical(
                f"Could not find analogue column in already existing scaffold.pkl at {reactant_analogues_path}. "
                f"Stopping...")

    def save(self):
        """
        This function saves the library object in the self.extra_output_dir_path.
        """
        with open(f"{self.extra_dir_path}/{self.id}_library.pkl", "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(output_dir: str, product: str) -> 'Library':
        """
        This function loads the library object. Will have to load all library objects
        """
        logger = logging.getLogger(f"{__name__}.{Library.__name__}")
        # get all .pkl in output_dir.
        library_pkls: List[str] = glob.glob(f"{output_dir}/*library*.pkl")
        # get id from scaffold SMILES
        id = fairy.generate_inchi_ID(product, isomeric=False)
        with open(library_pkls[0], "rb") as f:
            library = pickle.load(f)
            if library.id == id and library.num_steps == library.current_step:  # return the final library
                logger.info(f"Loaded {library.id} final library.")
                return library
        # if not found, return None
        logger.warning(f"Could not load {id} final library.")
        return None
