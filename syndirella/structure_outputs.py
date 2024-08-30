#!/usr/bin/env python3
"""
syndirella.structure_outputs.py

This module contains the functions used to structure the pipeline outputs.
"""
import logging
from typing import *
import glob2
import pandas as pd
import os
from rdkit import Chem
import traceback

from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper
import syndirella.fairy as fairy

logger = logging.getLogger(__name__)


def add_route_info(reaction_names: List[str],
                   reactants: List[Tuple[str]],
                   route_uuid: str) -> Dict:
    """
    This function formats the full route into seperate columns for each reactant.
    """
    num_steps: int = len(reaction_names)
    route: Dict[str, str] = {'route_uuid': route_uuid}
    for step in range(num_steps):
        route[f'{step + 1}_reaction'] = reaction_names[step]
        route[f'{step + 1}_r1_smiles'] = reactants[step][0]
        try:
            route[f'{step + 1}_r2_smiles'] = reactants[step][1]
        except IndexError:
            pass
    return route


def add_outcome_info(slipper: Slipper) -> Dict[str, Any]:
    """
    This function adds placement information (None if not attempted).
    """
    num_placed: int | None = None
    num_successful: int | None = None
    to_hippo_path: str | None = None
    template: str | None = None
    total_num_unique_products: int | None = None
    total_num_products_enumstereo: int | None = None
    hits_names: List[str] = slipper.hits_names
    try:
        template = slipper.template
        num_placed = slipper.num_placed
        num_successful = slipper.num_successful
        to_hippo_path = slipper.to_hippo_path
        total_num_unique_products = slipper.num_unique_products
        total_num_products_enumstereo = slipper.num_products_enumstereo
    except AttributeError:
        pass
    outcome: Dict = {'total_num_unique_products': total_num_unique_products,
                     'total_num_products_enumstereo': total_num_products_enumstereo,
                     'num_placed': num_placed,
                     'num_successful': num_successful,
                     'to_hippo': to_hippo_path,
                     'template': template}
    for i, hit in enumerate(hits_names):
        outcome[f"hit{i + 1}"] = hit
    return outcome


def get_output_df(csv_path: str,
                  output_dir: str) -> Tuple[pd.DataFrame, str | None]:
    """
    Given the csv_path and output_dir, checks and reads in the most recent previous output df.

    Format of output csv name: [name_of_input_csv]_output_YYYYMMDD_HHMM.csv
    """
    past_csv_path: str | None = None
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Could not find csv file at {csv_path}")
    csv_name: str = os.path.basename(csv_path).split('.')[0]
    # look for a csv that contains csv name in output dir
    output_csvs: List[str] = glob2.glob(os.path.join(output_dir, f'*{csv_name}_output*'))
    if len(output_csvs) == 0:
        # does not exist yet, make new blank df
        df = pd.DataFrame(columns=['smiles', 'inchi_key', 'route_uuid', 'error_type', 'error_message',
                                   'num_placed', 'num_successful', '1_reaction', '1_r1_smiles', 'hit1', 'template',
                                   'to_hippo'])
        return df, past_csv_path
    else:
        # Sort files by their date and time in the filename
        output_csvs.sort(
            key=lambda x: os.path.basename(x).split('_')[-2] + os.path.basename(x).split('_')[-1].split('.')[0],
            reverse=True)
    # make sure pandas can read in the csv
    try:
        df = pd.read_csv(output_csvs[0])
        past_csv_path: str = output_csvs[0]
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=['smiles', 'inchi_key', 'route_uuid', 'error_type', 'error_message',
                                   'num_placed', 'num_successful', '1_reaction', '1_r1_smiles', 'hit1', 'template',
                                   'to_hippo'])
    return df, past_csv_path


def get_scaffold_smiles(error: Exception | None, smiles: str | None, workshop: CobblersWorkshop | None) -> str:
    """
    This function gets the scaffold smiles from the error, smiles, or workshop.
    """
    if error is not None:  # highest level of detail
        try:
            scaffold = Chem.MolToSmiles(error.mol)
            return scaffold
        except AttributeError:
            try:
                scaffold = error.smiles
                return scaffold
            except AttributeError:
                pass
    if workshop is not None:  # next level of detail
        try:
            scaffold = workshop.product
            return scaffold
        except AttributeError:
            pass
    if smiles is not None:  # lowest level of detail
        return smiles
    else:
        raise ValueError("No scaffold found.")


def get_error_info(error: Exception | None) -> Tuple[str | None, str | None, bool]:
    """
    This function gets the error type and message from the error.
    """
    if error is not None:
        error_type: str | None = type(error).__name__
        try:
            error_message: str | None = error.message  # custom error
            custom_error: bool = True
        except AttributeError:
            error_message = error.args[0]  # any other error
            custom_error = False
    else:
        error_type = None
        error_message = None
        custom_error = False
    return error_type, error_message, custom_error


def get_inchi(scaffold: str, workshop: CobblersWorkshop | None) -> str:
    """
    This function gets the inchi from the scaffold and workshop.
    """
    if workshop is not None:
        try:
            return workshop.id
        except AttributeError:
            pass
    else:
        return fairy.generate_inchi_ID(scaffold)


def check_route_to_add(workshop: CobblersWorkshop | None) -> bool:
    """
    This function checks if there is a route defined in the workshop object.

    Returns True if route is defined.
    Returns False if route is not defined.
    """
    if workshop is not None:
        try:
            route_uuid: str = workshop.route_uuid
            reaction_names: List[str] = workshop.reaction_names
            reactants: List[Tuple[str, str]] = workshop.reactants
            return True
        except AttributeError:
            return False
    else:
        return False


def check_placement_to_add(slipper: Slipper | None) -> bool:
    """
    This function checks if there is placement information in the slipper object.

    Returns True if placement information exists.
    Returns False if placement information does not exist.
    """
    if slipper is not None:
        try:
            hits: List[str] = slipper.hits_names  # just need to see if slipper contains the least amount of info
            return True
        except AttributeError:
            return False
    else:
        return False


def check_additional_info_to_add(slipper: Slipper | None) -> bool:
    """
    This function checks if there is additional information in the slipper object.
    """
    if slipper is not None:
        try:
            additional_info: Dict = slipper.additional_info
            return True
        except AttributeError:
            return False
    else:
        return False

def add_new_route_to_output_df(output_df: pd.DataFrame, row: Dict) -> pd.DataFrame:
    """
    This function adds a new row to the output dataframe. If the row contains different columns,
    they are added to the output_df, and all values are filled with None if not present before.
    """
    # Convert row dictionary to DataFrame for easier manipulation
    row_df = pd.DataFrame([row])
    # Align the columns of output_df and row_df
    combined_columns = output_df.columns.union(row_df.columns)
    output_df = output_df.reindex(columns=combined_columns)
    row_df = row_df.reindex(columns=combined_columns)
    # Check for duplicate entries
    for i, r in output_df.iterrows():
        if r['smiles'] == row['smiles'] and r['route_uuid'] == row['route_uuid']:
            logger.info(f"Route already exists in output csv, replacing.")
            output_df.drop(i, inplace=True)
    # Append the new row to output_df
    output_df = pd.concat([output_df, row_df], ignore_index=True)

    # order columns with smiles, inchi_key, error_type, error_message always first
    try:
        output_df = output_df[['smiles', 'inchi_key', 'error_type', 'error_message'] + [col for col in output_df.columns if col not in ['smiles', 'inchi_key', 'error_type', 'error_message']]]
    except KeyError:
        pass
    return output_df

def save_output_df(output_df: pd.DataFrame, output_dir: str, csv_path: str):
    """
    This function saves the output dataframe to a csv file in the output directory.

    Format of output csv name: [name_of_input_csv]_output_YYYYMMDD_HHMM.csv
    """
    csv_name: str = os.path.basename(csv_path).split('.')[0]
    output_name: str = f'{csv_name}_output_{pd.Timestamp.now().strftime("%Y%m%d_%H%M")}.csv'
    output_csv_path: str = os.path.join(output_dir, output_name)
    output_df.to_csv(output_csv_path, index=False)
    logger.info(f"Output csv saved to {output_csv_path}")


def structure_route_outputs(error_message: str | None, error_type: str | None, output_df: pd.DataFrame,
                            workshop: CobblersWorkshop | None, slipper: Slipper | None, scaffold: str) -> pd.DataFrame:
    """
    This function structures the route outputs as a single row in the output dataframe.
    """
    inchi: str = get_inchi(scaffold=scaffold, workshop=workshop)
    row: Dict = {'smiles': scaffold,
                 'inchi_key': inchi,
                 'error_type': error_type,
                 'error_message': error_message}
    if check_route_to_add(workshop):  # True if route exists
        reaction_info: Dict | None = add_route_info(reaction_names=workshop.reaction_names,
                                                    reactants=workshop.reactants,
                                                    route_uuid=workshop.route_uuid)
        row.update(reaction_info)
    if check_placement_to_add(slipper):  # True if slipper exists and at least contains hits_names
        placement_info: Dict = add_outcome_info(slipper)
        row.update(placement_info)
    if check_additional_info_to_add(slipper):  # True if slipper contains additional info
        additional_info: Dict = slipper.additional_info
        row.update(additional_info)
    output_df: pd.DataFrame = add_new_route_to_output_df(output_df=output_df, row=row)
    return output_df


##########################################################

def structure_pipeline_outputs(error: Exception | None,
                               csv_path: str,
                               output_dir: str,
                               workshop: CobblersWorkshop | None = None,
                               slipper: Slipper | None = None,
                               smiles: str | None = None):
    """
    Structure outputs of pipeline.
    """
    try:
        output_df, past_csv_path = get_output_df(csv_path=csv_path, output_dir=output_dir)
        error_type, error_message, custom_error = get_error_info(error=error)
        if custom_error:
            scaffold: str = get_scaffold_smiles(error=error,
                                                smiles=smiles,
                                                workshop=workshop)
        else:
            scaffold: str = get_scaffold_smiles(error=None,
                                                smiles=smiles,
                                                workshop=workshop)
        output_df: pd.DataFrame = structure_route_outputs(error_message=error_message,
                                                          error_type=error_type,
                                                          output_df=output_df,
                                                          workshop=workshop,
                                                          slipper=slipper,
                                                          scaffold=scaffold)
        save_output_df(output_df=output_df, output_dir=output_dir, csv_path=csv_path)
        if past_csv_path is not None:
            os.remove(past_csv_path) # delete previous output csv
            logger.info(f"Deleted previous output csv at {past_csv_path}")
    except (TypeError, ValueError, FileNotFoundError):
        logger.error(f"Could not structure pipeline outputs.")
        logger.error(traceback.format_exc())
        # don't raise error, just log it


"""
TODO:
    Handle different structures of product (mol, smiles, inchi). 
    Handle message to put in error output. 
    Handle different levels of route matching:
        1. route_uuid
        2. inchi
        3. mol
        4. smiles
    Handle any different error. 
        Can still get message from error if not custom. 

Each row is a route. 

Cases to handle:
1. First finished route --> make new output csv
2. Not first finished route --> read in previous output csv, save new one

Columns (required)
    smiles: 
        Error.mol
        Error.smiles
        smiles (if not None)
        --> cannot get smiles from inchi
    inchi_key: 
        
    route_uuid: 
    error_type: None if hippo df is saved
    error_message: None if hippo df is saved
    num_placed: None if not attempted
    num_successful: None if not attempted / failed halfway
        Defined as ddG <0 , RMSD <2, passed PB
    reaction_name_step1: 
    reactant_step1:
    reactant2_step1:
    hit1:
    template: path
    [additional_info]:

# last step ends at reactant

Columns (not required):
    reaction_name_stepX:
    reactant_stepX:
    product_stepX:
    hitX: 

Places when this will be created:
    Fully successful, hippo df saved
    Failed in placement
    Failed by finding no stereoisomers
    Failed by finding no analogues --> no products
    Failed by SMARTS
    Failed by reaction
    Failed by creating route
    Failed by finding no routes
    Failed by couldn't read in smiles to mol
"""
