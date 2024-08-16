#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""
import os
from typing import List, Tuple
import pandas as pd
from rdkit import Chem
import datetime, time
import traceback
import glob2
import logging
from rdkit.Chem import inchi

from syndirella.Cobbler import Cobbler
from syndirella.cobblers_workshop.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper
from syndirella.slipper.SlipperFitter import SlipperFitter
import syndirella.fairy as fairy
from syndirella.error import ScaffoldPlacementError, NoSynthesisRoute, NoReactants
import syndirella.check_inputs as check_inputs

logger = logging.getLogger(__name__)

# TODO: Check if these functions are still needed after making check_inputs.py

def get_template_path(template_dir: str, template: str) -> str:
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


def assert_scaffold_placement(scaffold: str,
                               template_path: str,
                               hits_path: str,
                               hits_names: List[str],
                               output_dir: str
                               ) -> None:
    """
    Assert that the scaffold can be placed.
    """
    scaffold_mol = Chem.MolFromSmiles(scaffold)
    if not scaffold_mol:
        logger.critical(f"Could not create a molecule from the smiles {scaffold}.")
        raise ValueError(f"Could not create a molecule from the smiles {scaffold}.")
    slipper_fitter = SlipperFitter(template_path, hits_path, hits_names, output_dir)
    can_be_placed: bool = slipper_fitter.check_scaffold(scaffold_mol)
    if not can_be_placed:
        logger.critical(f"Scaffold {scaffold} could not be placed successfully.")
        raise ScaffoldPlacementError(scaffold)


def elaborate_from_cobbler_workshops(cobbler_workshops: List[CobblersWorkshop],
                                      template_path: str,
                                      hits_path: str,
                                      hits: List[str],
                                      batch_num: int,
                                      additional_info=None):
    """
    Does elaboration once the cobbler workshops are created.
    """
    if additional_info is None:
        additional_info = []
    for workshop in cobbler_workshops:
        try:
            final_library = workshop.get_final_library()
            if final_library is None:
                logger.warning(f"Could not get the final library for compound {workshop.product}. Skipping...")
                continue
            slipper = Slipper(library=final_library, template=template_path, hits_path=hits_path, hits_names=hits,
                              batch_num=batch_num, additional_info=additional_info)
            _, uuid = slipper.get_products()
            slipper.place_products()
            slipper.write_products_to_hippo(uuid=uuid)  # only write at the end after placement, to get correct route_uuid
            slipper.clean_up_placements()
        except Exception as e:
            tb = traceback.format_exc()
            logger.error(f"Error elaborating compound {workshop.product}. {tb}")
            continue


# def get_alternative_routes(product: str,
#                            reaction_names: List[str],
#                            reactants: List[Tuple[str, str]],
#                            num_steps: int,
#                            output_dir: str) -> CobblersWorkshop | None:
#     if fairy.do_i_need_alterative_route(reaction_names):
#         logger.info(f"Found the need for an alternative route for compound {product}.")
#         try:
#             cobbler = Cobbler(scaffold_compound=product, output_dir=output_dir)
#             cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
#             cobbler_workshops = [workshop for workshop in cobbler_workshops if
#                                  workshop.reaction_names != reaction_names]
#             workshop = CobblersWorkshop(product=product, reactants=reactants, reaction_names=reaction_names,
#                                         num_steps=num_steps, output_dir=output_dir, filter=True)
#         except Exception as e:
#             tb = traceback.format_exc()
#             logger.error(f"Error finding alternative route for compound {product}. {tb}")

def start_elaboration(product: str,
                      template_path: str,
                      hits_path: str,
                      hits: List[str],
                      output_dir: str,
                      additional_info=None):
    """
    Starts the elaboration of a single compound by checking if it can be placed.
    """
    if additional_info is None:
        additional_info = []
    logger.info(f'Running pipeline for: {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))}')
    start_time = time.time()
    mol = Chem.MolFromSmiles(product)
    if not mol:
        logger.error(f"Could not create a molecule from the smiles {product}.")
        raise ValueError(f"Could not create a molecule from the smiles {product}.")
    assert_scaffold_placement(scaffold=product, template_path=template_path, hits_path=hits_path, hits_names=hits,
                              output_dir=output_dir)
    return start_time

def elaborate_compound_with_manual_routes(product: str,
                                           reactants: List[Tuple[str, str]],
                                           reaction_names: List[str],
                                           num_steps: int,
                                           hits: List[str],
                                           template_path: str,
                                           hits_path: str,
                                           batch_num: int,
                                           output_dir: str,
                                           additional_info=None):
    """
    This function is used to elaborate a single compound using a manually defined route.
    """
    start_time: float = start_elaboration(product=product, template_path=template_path, hits_path=hits_path,
                                          hits=hits, output_dir=output_dir, additional_info=additional_info)
    workshop = CobblersWorkshop(product=product, reactants=reactants, reaction_names=reaction_names,
                                num_steps=num_steps, output_dir=output_dir, filter=False)
    cobblers_workshops = [workshop]
    alternative_routes: List[CobblersWorkshop] | None = workshop.get_additional_routes(edit_route=True)
    if alternative_routes is not None:
        cobblers_workshops = [workshop] + alternative_routes
    elaborate_from_cobbler_workshops(cobbler_workshops=cobblers_workshops, template_path=template_path,
                                     hits_path=hits_path, hits=hits, batch_num=batch_num,
                                     additional_info=additional_info)
    end_time = time.time()
    elapsed_time = end_time - start_time
    logger.info(f"Finished elaborating compound {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))} "
                f"after {datetime.timedelta(seconds=elapsed_time)}")
    logger.info("")


def elaborate_compound_full_auto(product: str,
                                  hits: List[str],
                                  template_path: str,
                                  hits_path: str,
                                  batch_num: int,
                                  output_dir: str,
                                  additional_info=None):
    """
    This function is used to elaborate a single compound.
    """
    start_time: float = start_elaboration(product=product, template_path=template_path, hits_path=hits_path,
                                          hits=hits, output_dir=output_dir, additional_info=additional_info)
    cobbler = Cobbler(scaffold_compound=product, output_dir=output_dir) # check that output_dirs can be made for different routes
    cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
    elaborate_from_cobbler_workshops(cobbler_workshops=cobbler_workshops, template_path=template_path,
                                     hits_path=hits_path, hits=hits, batch_num=batch_num,
                                     additional_info=additional_info)
    end_time = time.time()
    elapsed_time = end_time - start_time
    logger.info(f"Finished elaborating compound {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))} "
                f"after {datetime.timedelta(seconds=elapsed_time)}")
    logger.info("")


def run_pipeline(*,
                 csv_path: str,
                 output_dir: str,
                 template_dir: str,
                 hits_path: str,
                 metadata_path: str,
                 batch_num: int,
                 additional_columns: List[str],
                 manual_routes: bool,
                 ):
    """
    Run the whole syndirella pipeline! 👑
    """
    check_inputs.check_pipeline_inputs(csv_path=csv_path, template_dir=template_dir, hits_path=hits_path,
                                       metadata_path=metadata_path, additional_columns=additional_columns,
                                       manual_routes=manual_routes)
    df = pd.read_csv(csv_path)

    if not manual_routes:
        logger.info("Running the full auto pipeline.")
        for index, row in df.iterrows():
            try:
                additional_info: dict = check_inputs.format_additional_info(row, additional_columns)
                template_path: str = check_inputs.get_template_path(template_dir=template_dir, template=row['template'],
                                                                    metadata_path=metadata_path)
                hits: List[str] = check_inputs.get_exact_hit_names(row=row, metadata_path=metadata_path)
                elaborate_compound_full_auto(product=row['smiles'], hits=hits, template_path=template_path,
                                             hits_path=hits_path, batch_num=batch_num, output_dir=output_dir,
                                             additional_info=additional_info)
            except ScaffoldPlacementError:
                logger.warning(f"scaffold compound {row['smiles']} could not be placed successfully. Skipping...")
                continue
            except NoSynthesisRoute:
                logger.warning(f"No synthesis route found for compound {row['smiles']}. Skipping...")
                continue
            except NoReactants:
                logger.warning(f"No reactants found for compound {row['smiles']}. Skipping...")
                continue
    else:
        logger.info("Running the pipeline with manual routes.")
        for index, row in df.iterrows():
            try:
                additional_info: dict = check_inputs.format_additional_info(row, additional_columns)
                template_path: str = check_inputs.get_template_path(template_dir=template_dir, template=row['template'],
                                                                    metadata_path=metadata_path)
                reactants, reaction_names, num_steps = check_inputs.format_manual_route(row)
                hits: List[str] = check_inputs.get_exact_hit_names(row=row, metadata_path=metadata_path)
                elaborate_compound_with_manual_routes(product=row['smiles'], reactants=reactants,
                                                      reaction_names=reaction_names, num_steps=num_steps, hits=hits,
                                                      template_path=template_path, hits_path=hits_path,
                                                      batch_num=batch_num, output_dir=output_dir,
                                                      additional_info=additional_info)
            except ScaffoldPlacementError:
                logger.warning(f"scaffold compound {row['smiles']} could not be placed successfully. Skipping...")
                continue
            except NoSynthesisRoute:
                logger.warning(f"No synthesis route found for compound {row['smiles']}. Skipping...")
                continue
            except NoReactants:
                logger.warning(f"No reactants found for compound {row['smiles']}. Skipping...")
                continue

    logger.info("Pipeline complete.")
