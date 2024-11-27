#!/usr/bin/env python3
"""
syndirella.fairy.py

This module provides functions to output retrosynthesis queries for a given list of scaffolds.
"""
import os.path
from typing import Dict, List
import pandas as pd
import logging
from .cli_defaults import cli_default_settings
from syndirella.Postera import Postera
import json

logger = logging.getLogger(__name__)
with open(cli_default_settings['rxn_smarts_path']) as f:
    reaction_smarts = json.load(f)
reaction_smarts_names: List[str] = list(reaction_smarts.keys())

def save_df(df: pd.DataFrame, output_dir: str, csv_path: str) -> str:
    """
    Save the DataFrame to the output directory.
    """
    csv_basename = os.path.basename(csv_path)
    pkl_basename = csv_basename.replace('.csv', '.pkl.gz')
    saved_path = os.path.join(output_dir, f'justretroquery_{pkl_basename}')
    df.to_pickle(saved_path)
    return saved_path

def format_routes(routes: List[Dict[str, List[Dict[str, str]]]]) -> Dict:
    """
    Gets the top 5 passing routes from the routes. Formats them into a dictionary with routes names as keys, also adds
    other field of routeX_names.
    """
    passing_routes = {}
    n_of_rxns: List[int] = [i for i, route in enumerate(routes) if len(route['reactions']) > 0]
    if len(n_of_rxns) == 0:
        logger.error("No routes with reactions found.")
        return passing_routes
    min_i = min(n_of_rxns)
    for i, route in enumerate(routes):
        # get the reactions
        reactions: List[Dict[str, str]] = route['reactions']
        # check if all reactions are in the reactions
        if len(reactions) == 0:
            continue
        # adjust i to be in range(5)
        adjusted_i = i - min_i
        # in original reactions list, replace spaces with underscores
        for reaction in reactions:
            reaction['name'] = reaction['name'].replace(" ", "_")
        passing_routes[f'route{adjusted_i}'] = reactions
        passing_routes[f'route{adjusted_i}_names']: List[str] = [reaction['name'] for reaction in reactions]
        passing_routes[f'route{adjusted_i}_CAR']: bool = all(
            [True if reaction['name'] in reaction_smarts_names else False for reaction in reactions])
        passing_routes[f'route{adjusted_i}_non_CAR']: List[str] | None = ([reaction['name'] for reaction in reactions if
                                                                          reaction['name'] not in reaction_smarts_names]
                                                                          or None)
        logger.info(f"Route {adjusted_i} has {len(reactions)} reaction(s)")
        logger.info(f"Route {adjusted_i} CAR: {passing_routes[f'route{adjusted_i}_CAR']}")
        logger.info(f"Route {adjusted_i} non-CAR: {passing_routes[f'route{adjusted_i}_non_CAR']}")
        if len(passing_routes) == 20:  # 5 routes with 3 fields each
            break
    return passing_routes

def retro_search(scaffold: str) -> pd.DataFrame:
    """
    Perform retrosynthesis search on the given scaffold and formats outputs.
    """
    postera = Postera()
    routes = postera.perform_route_search(scaffold)
    formatted_routes: Dict = format_routes(routes)
    logger.info(len(formatted_routes))
    to_add = {'smiles': scaffold}
    for route, details in formatted_routes.items():
        to_add[route] = details
    return pd.DataFrame([to_add]) # each dictionary is a single row

def process_df(df: pd.DataFrame):
    """
    Process the input DataFrame and create output df with retrosynthesis information.
    """
    logger.info(f"Processing DataFrame with len {len(df)}")
    route_infos = []
    for i, scaffold in enumerate(df['smiles']):
        route_info: pd.DataFrame = retro_search(scaffold)
        route_infos.append(route_info)
    # format the DataFrame
    route_info_df = pd.concat(route_infos)
    merged_df = df.merge(route_info_df, on='smiles')
    merged_df.reset_index(drop=True, inplace=True)
    return merged_df

#######################################

def run_justretroquery(settings: Dict):
    """
    Run the justretroquery pipeline with the given settings.
    """
    # all you need is input csv with just a column of smiles
    try:
        csv_path: str = settings['input']
        output_dir: str = settings['output']
    except KeyError as e:
        raise KeyError(f"Missing critical argument to run justretroquery: {e}")

    df: pd.DataFrame = pd.read_csv(csv_path)

    df: pd.DataFrame = process_df(df)

    saved_path: str = save_df(df, output_dir, csv_path)

    logger.info(f"Saved DataFrame to {saved_path}")

    logger.info('Justretroquery execution completed successfully.')