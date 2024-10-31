#!/usr/bin/env python3
"""
syndirella.check_inputs.py

This module contains the functions used to check the inputs for running the pipeline.
"""
import os
import logging
import glob2
import pandas as pd
from typing import List, Dict, Set, Tuple, Any
import re

import numpy as np
from rdkit import Chem
from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)


def check_csv(csv_path: str) -> None:
    """
    Make sure that the csv path exists, can be read, and contains the required columns.
    """
    if not os.path.exists(csv_path):
        logger.error("The csv path does not exist.")
        raise FileNotFoundError(f"The csv path {csv_path} does not exist.")
    df = pd.read_csv(csv_path)
    required_columns = ['smiles', 'compound_set', 'template', 'hit1']  # must contain at least one hit
    for col in required_columns:
        if col not in df.columns:
            logger.critical(f"The csv must contain the column {col}.")
            raise ValueError(f"The csv must contain the column {col}.")


def metadata_dict(metadata_path: str, long_code_column: str = 'Long code') -> Dict:
    """
    Get the metadata dictionary from the metadata file, checking that it contains the required columns.
    """
    if not os.path.exists(metadata_path):
        logger.error("The metadata path does not exist.")
        raise FileNotFoundError(f"The metadata path {metadata_path} does not exist.")
    metadata = pd.read_csv(metadata_path)
    if 'Code' in metadata.columns and long_code_column in metadata.columns:
        code_info: Dict = metadata.set_index('Code')[long_code_column].to_dict()
    else:
        logger.warning(f"The metadata does not contain the columns 'Code' and {long_code_column}. Searching for 'crystal_name' "
                       "column instead.")
        if 'crystal_name' in metadata.columns:
            # use crystal_name as key and value to match dict
            logger.info("Using 'crystal_name' column as the key and value for the metadata dictionary.")
            code_info: Dict = {name: name for name in metadata['crystal_name']}
        else:
            raise ValueError(f"The metadata must contain the columns 'Code' and '{long_code_column}'.")
    return code_info


def check_template_paths(template_dir: str, csv_path: str, metadata_path: str) -> Set[str]:
    """
    Get the exact template paths, checking that they exist in the template directory.
    """
    if not os.path.isdir(template_dir):
        logger.error("The template directory does not exist.")
        raise NotADirectoryError(f"The template directory {template_dir} does not exist.")
    df = pd.read_csv(csv_path)
    templates: List[str] = [template.strip() for template in df['template'].tolist()]  # remove whitespace
    code_dict: Dict = metadata_dict(metadata_path)
    # get exact code for hit from metadata
    exact_codes = [key for key in code_dict for template in templates if template.lower() in key.lower()]
    template_paths = []
    for code in exact_codes:
        template_path = glob2.glob(f"{template_dir}/**/*{code}*.pdb")
        if len(template_path) == 0:
            logger.error(f"The template {code} does not exist in the template directory.")
            raise FileNotFoundError(f"The template {code} does not exist in the template directory.")
        elif len(template_path) > 1:
            logger.error(f"Multiple templates found for {code}. Please ensure that the template name is unique in the "
                         f"input csv.")
            raise ValueError(f"Multiple templates found for {code}. Please ensure that the template name is unique in"
                             f"the input csv.")
        template_paths.append(template_path[0])
    return set(template_paths)


def fill_in_product(row: pd.Series, step: int) -> None:
    """
    Fill in the scaffold for the given step.
    """
    # TODO
    # Use SMARTSHandler to do this
    pass


def check_route(i: int, row: pd.Series) -> None:
    """
    Checks that the route is in the correct format and fills in missing products if needed.
    """
    # fill series with nan if missing
    row = row.fillna(value='None')
    # find what number of steps the route is
    try:
        steps = [int(list(col.split('_')[-1])[-1]) for col in row.index if 'reaction_name' in col]
    except ValueError:
        logger.critical(f"Error in the route at row {i}. Please check names of columns to match the exact format"
                        f"as the template.")
        raise ValueError(f"Error in the route at row {i}. Please check names of columns to match the exact format"
                         f"as the template.")
    found = False
    last_step = max(steps)
    while found is False:
        if row[f'reaction_name_step{last_step}'] != 'None':  # use last filled reaction_name_stepX to find last step
            found = True
        else:
            last_step -= 1
    # check there are no missing reactants, reaction_names, or products for beginning and internal steps
    for step in range(1, last_step + 1):
        if step == 1:
            if row[f'reactant_step{step}'] == 'nan' or row[f'reactant2_step{step}'] == 'nan':  # only first step has 2
                logger.critical(f"Missing reactant for step {step} in route {i}.")
                raise ValueError(f"Missing reactant for step {step} in route {i}.")
        if step != last_step:
            if row[f'product_step{step}'] == 'nan':  # only internal steps have scaffold
                fill_in_product(row, step)
                logger.critical(f"Missing scaffold for step {step} in route {i}.")
                raise ValueError(f"Missing scaffold for step {step} in route {i}.")
        # check reaction_name
        if row[f'reaction_name_step{step}'] == 'nan':
            logger.critical(f"Missing reaction name for step {step} in route {i}.")
            raise ValueError(f"Missing reaction name for step {step} in route {i}.")
        # check reactant
        if row[f'reactant_step{step}'] == 'nan':
            logger.critical(f"Missing reactant for step {step} in route {i}.")
            raise ValueError(f"Missing reactant for step {step} in route {i}.")


def check_manual(csv_path: str) -> None:
    """
    Check that the manual dataframe is in the correct format, otherwise raise errors.
    """
    df = pd.read_csv(csv_path)
    # check that the route contains column names for at least 1 step
    required = ['reaction_name_step1', 'reactant_step1', 'reactant2_step1']
    if not all([col in df.columns for col in required]):
        logger.critical(f"Manual route must at least contain the columns {required}.")
        raise ValueError(f"Manual route must at least contain the columns {required}.")
    for i, row in df.iterrows():
        check_route(i, row)


def check_hit_names(csv_path: str, hits_path: str, metadata_path: str, long_code_column: str) -> None:
    """
    Check that the hit names are found within SDF.
    """
    df = pd.read_csv(csv_path)
    hit_cols = [col for col in df.columns if re.match(r'^hit\d+$', col)]
    hit_names = df[hit_cols].values.flatten()
    # remove nan and strip whitespace
    hit_names = [str(name).strip() for name in hit_names if str(name) != 'nan']
    if not os.path.exists(hits_path):
        logger.critical("The hits_path path does not exist.")
        raise FileNotFoundError(f"The hits_path path {hits_path} does not exist")
    sdf = Chem.SDMolSupplier(hits_path)
    sdf_names = [mol.GetProp('_Name') for mol in sdf]
    if not os.path.exists(metadata_path):
        logger.critical("The metadata path does not exist.")
        raise FileNotFoundError(f"The metadata path {metadata_path} does not exist.")
    code_dict = metadata_dict(metadata_path, long_code_column=long_code_column)
    # get the LongCodes for the hit names
    hit_longcodes = []
    for name in hit_names:
        matches = [code_dict[key] for key in code_dict if name.lower() in key.lower()]
        if len(matches) > 1:
            logger.critical(f"Multiple matches found in {metadata_path} using 'Code' for '{name}': {matches}. Please "
                            f"update the hit name in the input csv to be more specific.")
            raise ValueError(f"Multiple matches found for '{name}': {matches}")
        if len(matches) == 0:
            # could be an exact name for LongCode
            if name in sdf_names:
                hit_longcodes.append(name)
                continue
            logger.critical(f"No matches found in {metadata_path} using 'Code' for '{name}'. Please update the hit "
                            f"name in the input csv.")
            raise ValueError(f"No matches found for '{name}'")
        else:
            hit_longcodes.append(matches[0])
    # check if hit_longcodes are in sdf_names, not matching exactly, can be a substring
    if not all([any([longcode in sdf_name for sdf_name in sdf_names]) for longcode in hit_longcodes]):
        logger.critical(f"Not all hit names found in the sdf file. You might need to re-download hits_path and metadata"
                        f" from Fragalysis.")
        raise ValueError(
            f"Not all hit names found in the sdf file. You might need to re-download hits_path and metadata"
            f" from Fragalysis. Or set the code_column argument to the correct column in the metadata (such as "
            f"'Experiment code' or 'Compound code'.")


def check_apo_template(template_path: str) -> None:
    """
    Check that the template is actually apo (containing no LIG).
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('template', template_path)
    except Exception as e:
        logger.critical(f"Error parsing the template {template_path}: {e}")
        raise ValueError(f"Error parsing the template {template_path}: {e}")
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'LIG':
                    logger.critical(f"The template {template_path} contains 'LIG'. Please use an apo template.")
                    raise ValueError(f"The template {template_path} contains 'LIG'. Please use an apo template.")


## Additional functions for pipeline

def format_additional_info(row: pd.Series,
                           additional_columns: List[str]) -> Dict[str, Any]:
    """
    This function is used to format the additional info from the dataframe into a dictionary.
    """
    additional_info = {}
    for col in additional_columns:
        additional_info[col] = row[col]
    return additional_info


def get_exact_hit_names(row: pd.Series, metadata_path: str, hits_path: str) -> List[str]:
    """
    Get the exact hit name to use for placement.
    """
    code_dict = metadata_dict(metadata_path)
    hit_cols = [col for col in row.index if re.match(r'^hit\d+$', col)]
    hit_names = row[hit_cols].values.flatten()
    hit_names = [name.strip() for name in hit_names if str(name) != 'nan']
    sdf = Chem.SDMolSupplier(hits_path)
    sdf_names = [mol.GetProp('_Name') for mol in sdf]
    hit_longcodes = []
    for name in hit_names:
        matches = [code_dict[key] for key in code_dict if name.lower() in key.lower()]
        if len(matches) > 1:
            logger.critical(f"Multiple matches found in {metadata_path} using 'Code' for '{name}': {matches}. Please "
                            f"update the hit name in the input csv to be more specific.")
            raise ValueError(f"Multiple matches found in {metadata_path} using 'Code' for '{name}': {matches}. Please "
                             f"update the hit name in the input csv to be more specific.")
        if len(matches) == 0:
            # could be an exact name for LongCode
            if name in sdf_names:
                hit_longcodes.append(name)
                continue
            logger.critical(
                f"No matches found in {metadata_path} using 'Code' or in the hits SDF provided for '{name}'."
                f" Please update the hit name in the input csv.")
            raise ValueError(
                f"No matches found in {metadata_path} using 'Code' or in the hits SDF provided for '{name}'."
                f" Please update the hit name in the input csv.")
        else:
            hit_longcodes.append(matches[0])
    return hit_longcodes


def get_template_path(template_dir: str, template: str, metadata_path: str) -> str:
    """
    Get the exact template path to use for placement.
    """
    code_dict = metadata_dict(metadata_path)
    exact_code = [key for key in code_dict if template.lower() in key.lower()]
    if len(exact_code) == 0:
        logger.critical(f"The template {template} does not exist in the metadata.")
        raise ValueError(f"The template {template} does not exist in the metadata.")
    template_path = glob2.glob(f"{template_dir}/**/*{exact_code[0]}*.pdb")
    if len(template_path) == 0:
        logger.critical(f"The template {exact_code[0]} does not exist in the template directory.")
        raise FileNotFoundError(f"The template {exact_code[0]} does not exist in the template directory.")
    elif len(template_path) > 1:
        logger.critical(f"Multiple templates found for {exact_code[0]}. Please ensure that the template name is unique "
                        f"in the input csv.")
        raise ValueError(f"Multiple templates found for {exact_code[0]}. Please ensure that the template name is unique"
                         f"in the input csv.")
    return template_path[0]


def check_additional_columns(csv_path: str, additional_columns: List[str]) -> None:
    """
    Check that the additional columns are in the csv.
    """
    df = pd.read_csv(csv_path)
    for col in additional_columns:
        if col not in df.columns:
            logger.critical(f"The csv must contain the column {col} to add your desired metadata.")
            raise ValueError(f"The csv must contain the column {col} to add your desired metadata.")


def format_manual_route(row: pd.Series) -> Tuple[List[Tuple[Any, Any]], List[Any], int]:
    """
    Format route to output reactants, reaction names, and products.
    """
    # fill nan values
    row = row.fillna(value='None')
    # get the number of steps
    steps = [int(list(col.split('_')[-1])[-1]) for col in row.index if 'reaction_name' in col]
    found = False
    last_step: int = max(steps)
    while found is False:
        if row[f'reaction_name_step{last_step}'] == 'None':  # use last filled reaction_name_stepX to find last step
            last_step -= 1
        else:
            found = True
    # format reactants
    reactants = []
    reaction_names = []
    for step in range(1, last_step + 1):
        if step == 1:
            reactants_step = (row[f'reactant_step{step}'], row[f'reactant2_step{step}'])
            reactants.append(reactants_step)
        else:
            # one reactant is scaffold of previous step
            reactants_step = (row[f'product_step{step - 1}'], row[f'reactant_step{step}'])
            reactants.append(reactants_step)
        reaction_names.append(row[f'reaction_name_step{step}'])
    if len(reactants) != last_step:
        logger.critical(f"The number of reactants found does not match the number of steps in route for {row}.")
        raise ValueError(f"The number of reactants found does not match the number of steps in route for {row}.")
    if len(reaction_names) != last_step:
        logger.critical(f"The number of reaction names found does not match the number of steps in route for {row}.")
        raise ValueError(f"The number of reaction names found does not match the number of steps in route for {row}.")
    return reactants, reaction_names, last_step


###############################################

def check_pipeline_inputs(*,
                          csv_path: str,
                          template_dir: str,
                          hits_path: str,
                          metadata_path: str,
                          additional_columns: List[str],
                          manual_routes: bool,
                          long_code_column: str) -> None:
    """
    Check the inputs for the pipeline.
    """
    try:
        check_csv(csv_path)
        check_hit_names(csv_path, hits_path, metadata_path, long_code_column=long_code_column)
        check_additional_columns(csv_path, additional_columns)
        template_paths: Set[str] = check_template_paths(template_dir, csv_path, metadata_path)
        for template_path in template_paths:  # check each template
            check_apo_template(template_path)
        if manual_routes:
            check_manual(csv_path)
    except (TypeError, ValueError, FileNotFoundError):
        raise
    logger.info(f'All inputs are valid for {csv_path}.')
