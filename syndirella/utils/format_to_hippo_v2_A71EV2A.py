import argparse
import re
import os
import csv
import json
import pandas as pd
import shutil
import traceback

"""
This script will find the succesful placements of Fragmenstein. All that is needed in each compound dir are:
INPUT:
1. home_dir/cmpd_catalog_frag1_frag2/cmpd_catalog_frag1_frag2.csv (holds catalogue metadata for elaborations)
2. home_dir/cmpd_catalog_frag1_frag2/output (holds the output files from Fragmenstein)
NOTE: Does not require output.csv from Fragmenstein since sometimes full placements are not always completed. 

OUTPUT:
1. home_dir/cmpd_catalog_frag1_frag2/success (holds the successful placements)
2. home_dir/cmpd_catalog_frag1_frag2/success/cmpd_catalog_frag1_frag2_success.csv (holds the metadata for the successful placements)
3. /home_dir/success_dirs.json 

"""

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', required=True, help='Home directory for all the placements. '
                                                                 'Can be home directory of batches.')
    parser.add_argument('--rmsd', required=False, help='RMSD threshold to accept placement. All placements '
                                                       'below this will be accepted. Default is 2.0.')
    parser.add_argument('--remove', action='store_true', required=False,
                        help='Remove the extra fragmenstein files when moving'
                             'to success directory. Default is False.')
    parser.add_argument('--num_steps', required=True, help='Number of steps in route to help with qualifier addition.')
    return parser

# -------------------------#

def add_placement_data(elab_csv_path, success_csv_path):
    # Read elab csv
    elab_df = pd.read_csv(elab_csv_path)
    elab_df['name'] = elab_df['name'].str.replace('_', '-')
    # Read success csv
    success_df = pd.read_csv(success_csv_path)
    # Merge elab_df with success_df on 'name' column
    merged_df = pd.merge(success_df, elab_df, on='name', how='left')
    # drop index
    if 'index' in merged_df.columns:
        merged_df = merged_df.drop(columns=['index'])
    # order merged_df by 'num_atom_difference' column
    merged_df = merged_df.sort_values(by=['num_atom_difference'])
    # Write to success csv
    merged_df.to_csv(success_csv_path, index=False)
    # Write elab_df to elab csv
    elab_df.to_csv(elab_csv_path, index=False)

def find_num_placement(output_dir_path):
    # Function to count the number of first-level subdirectories in output_dir_path
    num_subdirs = 0
    if os.path.exists(output_dir_path):
        for entry in os.listdir(output_dir_path):
            if os.path.isdir(os.path.join(output_dir_path, entry)):
                num_subdirs += 1
    return num_subdirs

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

def read_json_and_check_conditions(json_path, rmsd_thresh):
    with open(json_path, 'r') as file:
        data = json.load(file)
        mRMSD = data.get("mRMSD", float('inf'))
        delta_delta_G = get_delta_delta_G(data)
        if mRMSD < rmsd_thresh and delta_delta_G < 0:
            return data
    return None

def append_to_csv(file_path, row_data):
    with open(file_path, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(row_data)

def copy_successful_cases(output_dir_path, success_dir_path, rmsd_thresh, remove=False):
    # Ensure success directory exists
    if not os.path.exists(success_dir_path):
        os.makedirs(success_dir_path)
    # Copy successful cases to success_dir
    for subdir in os.listdir(output_dir_path):
        if subdir.startswith('.'):  # Skip hidden files/directories
            continue
        subdir_path = os.path.join(output_dir_path, subdir)
        if os.path.isdir(subdir_path):
            for file in os.listdir(subdir_path):
                if file.endswith('.json'):
                    json_file_path = os.path.join(subdir_path, file)
                    if os.path.exists(json_file_path):
                        data = read_json_and_check_conditions(json_file_path, rmsd_thresh)
                        if data:
                            if remove:
                                # remove extra files in subdir_path that don't have these suffixes
                                suffixes_to_keep = ['.minimised.json', '.holo_minimised.pdb', '.minimised.mol', '.csv']
                                for file in os.listdir(subdir_path):
                                    if not any(file.endswith(suffix) for suffix in suffixes_to_keep):
                                        file_path = os.path.join(subdir_path, file)
                                        try:
                                            if os.path.isfile(file_path) or os.path.islink(file_path):
                                                os.unlink(file_path)
                                                print(f"Deleted file: {file_path}")
                                        except Exception as e:
                                            print('Failed to delete %s. Reason: %s' % (file_path, e))
                            # Copy the directory, handle existing directory case
                            destination_path = os.path.join(success_dir_path, subdir)
                            try:
                                shutil.copytree(subdir_path, destination_path)
                            except FileExistsError:
                                # Remove existing directory and re-copy
                                shutil.rmtree(destination_path)
                                shutil.copytree(subdir_path, destination_path)
                            break


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

def create_success_csv(success_dir_path):
    headers = ['name', 'ΔΔG', 'unbound', 'bound', 'mRMSD']
    collected_data = []
    for subdir in os.listdir(success_dir_path):
        if subdir.startswith('.'):  # Skip hidden files/directories
            continue
        subdir_path = os.path.join(success_dir_path, subdir)
        if os.path.isdir(subdir_path):
            for file in os.listdir(subdir_path):
                if file.endswith('.json'):
                    json_file_path = os.path.join(subdir_path, file)
                    with open(json_file_path, 'r') as file:
                        data = json.load(file)
                        try:
                            bound_value = data["Energy"]['bound']
                            placement_method = 'rdkit'
                        except KeyError:
                            placement_method = 'pyrosetta'
                        csv_row = make_success_csv_row(subdir, data)
                        collected_data.append(csv_row)
    csv_file_path = os.path.join(success_dir_path, f'{success_dir_path.split("/")[-2]}_success.csv')
    if collected_data:
        num_success = len(collected_data)
        with open(csv_file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            writer.writerows(collected_data)
        return csv_file_path, placement_method, num_success
    return None, None, None

def output_dirs_match(output_path, cmpd_catalog, frag1, frag2):
    """ Check if the output directory names match the cmpd_catalog, frag1, frag2. """
    output_dirs = os.listdir(output_path)
    cmpd_catalog = cmpd_catalog.replace('_','-')
    frag1 = frag1.replace('_','-')
    frag2 = frag2.replace('_','-')
    for output_dir in output_dirs:
        if cmpd_catalog in output_dir and frag1 in output_dir and frag2 in output_dir:
            return True
    return False

def extract_info(directory):
    """ Extracts the cmpd_catalog, frag1, and frag2 from the directory name. """
    pattern = r'(.+?)_(x\d+_\d+[AB])_(x\d+_\d+[AB])'
    match = re.search(pattern, directory)
    if match:
        return match.groups()
    else:
        return None, None, None

def find_cmpd_dirs(home_directory):
    """ Find directories that match the specified pattern. """
    cmpd_dirs = []
    for root, directories, _ in os.walk(home_directory):
        for directory in directories:
            if directory.startswith('.'):  # Skip hidden directories
                continue
            if 'batch' not in directory and 'output' not in directory and 'success' not in directory and 'results' not in directory:
                full_dir_path = os.path.join(root, directory)
                cmpd_dirs.append(full_dir_path)
    return cmpd_dirs

def contains_elab_csv(directory: str, num_steps: int):
    """Check if directory contains a .csv file with specific names."""
    if num_steps == 1:
        suffix = '1_of_1'
    else:
        suffix = '2_of_2'
    for file in os.listdir(directory):
        if file.startswith('.'):  # Skip hidden files
            continue
        if (file.endswith(".csv") and 'batch' not in file and 'success' not in file and suffix in file):
            elab_file_path = os.path.join(directory, file)
            return elab_file_path, True
    return None, False

# ----------------------------#
def main():
    parser = config_parser()
    args = parser.parse_args()
    total_success_dirs = []
    cmpd_dirs = find_cmpd_dirs(args.directory)
    for directory in cmpd_dirs:
        if directory == 'success' or directory == 'output' or 'xtra_results' in directory:
            continue
        # directory = name of the elab compound
        try:
            print('DIRECTORY:', directory)
            if args.rmsd is None:
                rmsd_thresh = 2.0
            else: rmsd_thresh = float(args.rmsd)
            # 1. If output directory does not exist within dir, skip.
            output_dir_path = os.path.join(directory, 'output')
            if not os.path.exists(output_dir_path):
                continue
            # 3. If there is no cmpd_catalog, frag1, frag2 in .csv, skip.
            num_steps = int(args.num_steps)
            elab_csv_path, found = contains_elab_csv(directory, num_steps)
            if not found:
                print(f"No relevant .csv file found in {directory}")
                continue
            print(f"FOUND elab_csv_path: {elab_csv_path}")
            # 5. Make success dir if it does not already exist.
            success_dir_path = os.path.join(directory, 'success')
            if not os.path.exists(success_dir_path):
                os.makedirs(success_dir_path)
            # 6. Find + copy successful placements.
            copy_successful_cases(output_dir_path, success_dir_path, rmsd_thresh, args.remove)
            # 7. Save success.csv.
            success_csv_path, placement_method, num_success = create_success_csv(success_dir_path)
            if success_csv_path is None:
                print(f"{directory} DOES NOT CONTAIN SUCCESSFUL PLACEMENTS")
                continue
            # 7. Find total number of placements by getting num of subdirs in output_dir_path
            num_placements = find_num_placement(output_dir_path)
            # 8. Add placement data to elab.csv and success.csv
            add_placement_data(elab_csv_path, success_csv_path)
            # 9. Store number of successful placements and total number of placements per directory. Store success dir paths.
            total_success_dirs.append({success_dir_path: (num_success, num_placements)})
            print(f"SUCCESS: {directory}")
        except Exception as e:
            tb = traceback.format_exc()
            print(f"Error in {directory}:\n{tb}")

    # Save results
    with open((os.path.join(args.directory, 'success_dirs.json')), 'w') as f:
        json.dump(total_success_dirs, f)


if __name__ == "__main__":
    main()
