#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""
import os
from typing import List, Tuple, Dict, Any
import pandas as pd
from rdkit import Chem
import datetime, time
import ast
import traceback
import glob2
import logging

from syndirella.Cobbler import Cobbler
from syndirella.cobblers_workshop.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper
from syndirella.slipper.SlipperFitter import SlipperFitter
from syndirella.Fairy import Fairy
from syndirella.error import RouteError

logger = logging.getLogger(__name__)

def _format_additional_info(row: pd.Series,
                            additional_columns: List[str]) -> Dict[str, Any]:
    """
    This function is used to format the additional info from the dataframe into a dictionary.
    """
    additional_info = {}
    for col in additional_columns:
        additional_info[col] = row[col]
    return additional_info


def _assert_csv(csv_path: str) -> str:
    """
    Make sure that the csv path exists.
    """
    if not os.path.exists(csv_path):
        logger.error("The csv path does not exist.")
        raise FileNotFoundError(f"The csv path {csv_path} does not exist.")
    df = pd.read_csv(csv_path)
    required_columns = ['smiles', 'compound_set', 'hits', 'template']
    for col in required_columns:
        if col not in df.columns:
            logger.error(f"The csv must contain the column {col}.")
            raise ValueError(f"The csv must contain the column {col}.")
    for index, row in df.iterrows():
        hits: str = str(row['hits'])
        split_hits = hits.split(' ')
        if not all(isinstance(hit, str) for hit in split_hits):
            logger.error(f"The 'hits' column must contain space-separated strings at index {index}.")
            raise ValueError(f"The 'hits' column must contain space-separated strings at index {index}.")
    return csv_path


def _get_template_path(template_dir: str, template: str) -> str:
    """
    Get the template path from the template directory.
    """
    if not os.path.isdir(template_dir):
        logger.error("The template directory does not exist.")
        raise NotADirectoryError(f"The template directory {template_dir} does not exist.")
    template_path: list = glob2.glob(f'{template_dir}/**/*{template}*', recursive=True)
    if len(template_path) != 1:
        logger.error(f"Could not find the template {template} in the template directory.")
        raise FileNotFoundError(f"Could not find the template {template} in the template directory.")
    if not os.path.exists(template_path[0]):
        logger.error(f"The template {template} does not exist.")
        raise FileNotFoundError(f"The template {template} does not exist.")
    return template_path[0]


def _assert_manual_df(df: pd.DataFrame) -> None:
    """
    Assert that the manual dataframe is in the correct format.
    """
    required_columns = ['smiles', 'reactants', 'reaction_names', 'num_steps']
    for col in required_columns:
        if col not in df.columns:
            logger.error(f"If doing a manual route, the csv must contain the column {col}.")
            raise ValueError(f"If doing a manual route, the csv must contain the column {col}.")
    for index, row in df.iterrows():
        try:
            reactants = ast.literal_eval(row['reactants'])
            reaction_names = ast.literal_eval(row['reaction_names'])
            num_steps = row['num_steps']
        except (ValueError, SyntaxError) as e:
            logger.error(f"Error parsing columns at index {index}: {e}")
            raise ValueError(f"Error parsing columns at index {index}: {e}")
        if len(reactants) != len(reaction_names):
            logger.error("The reactants and reaction_names columns must have the same number of elements.")
            raise ValueError("The reactants and reaction_names columns must have the same number of elements.")
        if not isinstance(reactants, list):
            logger.error("The reactants column must be a list of strings.")
            raise TypeError("The reactants column must be a list of strings.")
        if not isinstance(reaction_names, list):
            logger.error("The reaction_names column must be a list of strings.")
            raise TypeError("The reaction_names column must be a list of strings.")
        if not isinstance(num_steps, int):
            logger.error("The num_steps column must be an integer.")
            raise TypeError("The num_steps column must be an integer.")


def _assert_base_placement(base: str,
                           template_path: str,
                           hits_path: str,
                           hits_names: List[str],
                           output_dir: str
                           ) -> None:
    """
    Assert that the base can be placed.
    """
    base_mol = Chem.MolFromSmiles(base)
    if not base_mol:
        logger.error(f"Could not create a molecule from the smiles {base}.")
        raise ValueError(f"Could not create a molecule from the smiles {base}.")
    slipper_fitter = SlipperFitter(template_path, hits_path, hits_names, output_dir)
    can_be_placed: bool = slipper_fitter.check_base(base_mol)
    if not can_be_placed:
        logger.error(f"Base {base} could not be placed successfully.")
        raise RouteError(base)


def _elaborate_from_cobbler_workshops(cobbler_workshops: List[CobblersWorkshop],
                                      template_path: str,
                                      hits_path: str,
                                      hits: str,
                                      batch_num: int,
                                      additional_info: Dict[str, Any] = []):
    """
    Does elaboration once the cobbler workshops are created.
    """
    for workshop in cobbler_workshops:
        try:
            final_library = workshop.get_final_library()
            if final_library is None:
                logger.warning(f"Could not get the final library for compound {workshop.product}. Skipping...")
                continue
            slipper = Slipper(library=final_library,
                              template=template_path,
                              hits=hits_path,
                              hits_names=hits,
                              batch_num=batch_num,
                              additional_info=additional_info)
            _, uuid = slipper.get_products()
            slipper.place_products()
            slipper.write_products_to_hippo(uuid=uuid)  # only write at the end after placement, to get correct uuid
            slipper.clean_up_placements()
        except Exception as e:
            tb = traceback.format_exc()
            logger.error(f"Error elaborating compound {workshop.product}. {tb}")
            continue


def _elaborate_compound_with_manual_routes(product: str,
                                           reactants: List[Tuple[str]],
                                           reaction_names: List[str],
                                           num_steps: int,
                                           hits: str,
                                           template_path: str,
                                           hits_path: str,
                                           batch_num: int,
                                           output_dir: str,
                                           additional_info: Dict[str, Any] = []):
    """
    This function is used to elaborate a single compound using a manually defined route.
    """
    start_time = time.time()
    fairy = Fairy()
    mol = Chem.MolFromSmiles(product)
    if not mol:
        logger.error(f"Could not create a molecule from the smiles {product}.")
        raise ValueError(f"Could not create a molecule from the smiles {product}.")
    logger.info(f'Elaborating compound: {product}')
    hits = hits.split()
    _assert_base_placement(base=product, template_path=template_path, hits_path=hits_path, hits_names=hits,
                           output_dir=output_dir)

    workshop = CobblersWorkshop(product=product, reactants=reactants, reaction_names=reaction_names,
                                num_steps=num_steps, output_dir=output_dir, filter=False)
    cobbler_workshops = []
    if fairy.do_i_need_alterative_route(reaction_names):
        logger.info(f"Found the need for an alternative route for compound {product}.")
        try:
            cobbler = Cobbler(base_compound=product, output_dir=output_dir)
            cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
            cobbler_workshops = [workshop for workshop in cobbler_workshops if
                                 workshop.reaction_names != reaction_names]
            workshop = CobblersWorkshop(product=product, reactants=reactants, reaction_names=reaction_names,
                                        num_steps=num_steps, output_dir=output_dir, filter=True)
        except Exception as e:
            tb = traceback.format_exc()
            logger.error(f"Error finding alternative route for compound {product}. {tb}")
            cobbler_workshops = []
    cobbler_workshops.insert(0, workshop)
    _elaborate_from_cobbler_workshops(cobbler_workshops=cobbler_workshops, template_path=template_path,
                                      hits_path=hits_path, hits=hits, batch_num=batch_num,
                                      additional_info=additional_info)
    end_time = time.time()
    elapsed_time = end_time - start_time
    logger.info(f"Finished elaborating compound {product} after {datetime.timedelta(seconds=elapsed_time)}")
    logger.info("")


def _elaborate_compound_full_auto(product: str,
                                  hits: str,
                                  template_path: str,
                                  hits_path: str,
                                  batch_num: int,
                                  output_dir: str,
                                  additional_info: Dict[str, Any] = []):
    """
    This function is used to elaborate a single compound.
    """
    start_time = time.time()
    mol = Chem.MolFromSmiles(product)
    if not mol:
        logger.error(f"Could not create a molecule from the smiles {product}.")
        raise ValueError(f"Could not create a molecule from the smiles {product}.")
    hits = hits.split()
    _assert_base_placement(base=product, template_path=template_path, hits_path=hits_path, hits_names=hits,
                           output_dir=output_dir)

    cobbler = Cobbler(base_compound=product, output_dir=output_dir)
    cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
    _elaborate_from_cobbler_workshops(cobbler_workshops=cobbler_workshops, template_path=template_path,
                                      hits_path=hits_path, hits=hits, batch_num=batch_num,
                                      additional_info=additional_info)
    end_time = time.time()
    elapsed_time = end_time - start_time
    logger.info(f"Finished elaborating compound {product} after {datetime.timedelta(seconds=elapsed_time)}")
    logger.info("")


def run_pipeline(*,
                 csv_path: str,
                 output_dir: str,
                 template_dir: str,
                 hits_path: str,
                 batch_num: int,
                 additional_columns: List[str],
                 manual_routes: bool,
                 ):
    """
    Run the whole syndirella pipeline! 👑
    """
    csv_path = _assert_csv(csv_path)
    df = pd.read_csv(csv_path)

    for col in additional_columns:
        if col not in df.columns:
            logger.error(f"The csv must contain the column {col}.")
            raise ValueError(f"The csv must contain the column {col}.")

    if not manual_routes:
        logger.info("Running the full auto pipeline.")
        for index, row in df.iterrows():
            try:
                additional_info = _format_additional_info(row, additional_columns)
                template_path: str = _get_template_path(template_dir=template_dir, template=row['template'])
                _elaborate_compound_full_auto(product=row['smiles'], hits=row['hits'], template_path=template_path,
                                              hits_path=hits_path, batch_num=batch_num, output_dir=output_dir,
                                              additional_info=additional_info)
            except RouteError:
                logger.warning(f"Base compound {row['smiles']} could not be placed successfully. Skipping...")
                continue
    else:
        logger.info("Running the pipeline with manual routes.")
        _assert_manual_df(df)
        for index, row in df.iterrows():
            try:
                additional_info = _format_additional_info(row, additional_columns)
                template_path: str = _get_template_path(template_dir=template_dir, template=row['template'])
                _elaborate_compound_with_manual_routes(product=row['smiles'],
                                                       reactants=ast.literal_eval(row['reactants']),
                                                       reaction_names=ast.literal_eval(row['reaction_names']),
                                                       num_steps=row['num_steps'], hits=row['hits'],
                                                       template_path=template_path, hits_path=hits_path,
                                                       batch_num=batch_num, output_dir=output_dir,
                                                       additional_info=additional_info)
            except RouteError:
                logger.warning(f"Base compound {row['smiles']} could not be placed successfully. Skipping...")
                continue

    logger.info("Pipeline complete.")
