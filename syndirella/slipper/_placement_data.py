#!venv/bin/env python3
"""
syndirella.slipper._placement_data.py

This module contains functions to get the placement data from a Fragmenstein run.
"""
import os
import json
import pandas as pd
import csv
import glob2

def get_placement_data(products: pd.DataFrame, fragmenstein_output: str, library_output_dir: str) -> pd.DataFrame:
    # Make fragmenstein_placements.csv
    placements_path: str = make_fragmenstein_placements_csv(fragmenstein_output, library_output_dir)
    # Read csv
    placements: pd.DataFrame = pd.read_csv(placements_path)
    # Merge products with success_df on 'name' column
    merged_df = pd.merge(products, placements, on='name', how='left')
    # drop index
    if 'index' in merged_df.columns:
        merged_df = merged_df.drop(columns=['index'])
    # order merged_df by 'num_atom_difference' column
    merged_df = merged_df.sort_values(by=['num_atom_diff'])
    if 'Unnamed: 0' in merged_df.columns:
        merged_df = merged_df.drop(columns=['Unnamed: 0'])
    final_products_csv_path = find_products_csv(library_output_dir)
    merged_output_path = final_products_csv_path.split('.csv')[0] + '_placements.csv'
    merged_df.to_csv(merged_output_path, index=False)
    return merged_df

def find_products_csv(library_output_dir: str) -> str:
    """
    Given a directory, find the products and non-placements csv file.
    """
    # Pattern for matching files that contain 'products' but not '_placements' in their names
    pattern = library_output_dir + '/*products*.csv'
    matched_files = glob2.glob(pattern)
    # Filtering out files that contain '_placements'
    filtered_files = [file for file in matched_files if '_placements' not in file]
    return filtered_files[0]

def make_fragmenstein_placements_csv(output_path: str, library_output_dir: str) -> str:
    """
    This function makes a fragmenstein_placements.csv by looking through outputs of Fragmenstein.

    Args:
        output_path: str: The path to the output directory of Fragmenstein.
        library_output_dir: str: The path to the directory where the final products csv is located.
    """
    # Make fragmenstein_placements.csv
    headers = ['name', 'ΔΔG', 'ΔG_bound', 'ΔG_unbound', 'comRMSD']
    collected_data = []
    for subdir in os.listdir(output_path):
        if subdir.startswith('.'):  # Skip hidden files/directories
            continue
        subdir_path = os.path.join(output_path, subdir)
        if os.path.isdir(subdir_path):
            for file in os.listdir(subdir_path):
                if file.endswith('.json'):
                    json_file_path = os.path.join(subdir_path, file)
                    with open(json_file_path, 'r') as file:
                        data = json.load(file)
                        try:
                            csv_row = make_success_csv_row(subdir, data)
                            collected_data.append(csv_row)
                        except Exception as e:
                            print(f'Error: {e}')
    csv_file_path = os.path.join(library_output_dir, 'fragmenstein_placements.csv')
    if collected_data:
        with open(csv_file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            writer.writerows(collected_data)
        return csv_file_path
    return None

def get_delta_delta_G(data: dict) -> float:
    # Get the delta delta G value from the JSON file. Accounts for different formats.
    try:
        return data["Energy"]["xyz_∆∆G"]
    except KeyError:
        try:
            bound = data["Energy"]["bound"]['total_score']
            unbound = data["Energy"]["unbound"]['total_score']
            ddG = bound - unbound
            return ddG
        except KeyError:
            return float('inf')

def get_bound_unbound(data: dict) -> tuple:
    """Get the bound and unbound energy values from the JSON file. Accounts for different formats"""
    try:
        bound = data["Energy"]["xyz_bound"]
        unbound = data["Energy"]["xyz_unbound"]
        return bound, unbound
    except KeyError:
        try:
            bound = data["Energy"]["bound"]["total_score"]
            unbound = data["Energy"]["unbound"]["total_score"]
            return bound, unbound
        except KeyError:
            return float('inf'), float('inf')

def make_success_csv_row(subdir: str, data: dict) -> list:
    """Make a row for the success.csv file."""
    ddG = get_delta_delta_G(data)
    bound, unbound = get_bound_unbound(data)
    try:
        rmsd = data["mRMSD"]
    except KeyError:
        rmsd = float('inf')
    csv_row = [subdir, ddG, unbound, bound, rmsd]
    return csv_row
