import argparse
import os
import csv
import json
from typing import Dict, Any, Tuple
import glob2
import pandas as pd
import shutil
import unicodedata
from pathlib import Path

def config_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', required=True, help='Home directory for all the placements')
    parser.add_argument('-e', '--elab', required=True, help='Identifier of the elab CSV file to merge')
    parser.add_argument('-o', '--output', required=True, help='Identifier of the output CSV file to merge')
    parser.add_argument('--rmsd', required=False, help='RMSD threshold to accept placement. All placements '
                                                       'below this will be accepted. Default is 2.0.')
    parser.add_argument('--remove', action='store_true', required=False,
                        help='Remove the extra fragmenstein files when moving'
                             'to success directory. Default is False.')
    parser.add_argument('--orig_elab', required=True, help='Identifier of where to look for original elab csv. '
                                                           'Done if elab csv does not exist if it was accidentally deleted.')
    return parser


# -------------------------#

def find_create_elab_csv(root_path: str, dir_path: str, orig_elabs: str):
    """
    Creates elab.csv if it does not exist and there is an output folder. Should name it differently to avoid
    overwriting the original elab.csv.

    INPUT:
        root_path: root path to parent directory
        dir_path: directory name
    """
    print('elab_identifier', elab_identifier)
    print('output_identifier', output_identifier)
    elab_path = ''
    output_path = ''
    full_path = os.path.join(root_path, dir_path)
    print('full_path', full_path)
    # Find all csv files in the directory
    csv_files = glob2.glob(os.path.join(full_path, '*.csv'))
    print('csv_files', csv_files)
    # Find the elab and output files based on the unique identifiers
    for file in csv_files:
        if elab_identifier in file and output_identifier not in file and 'success_moved' not in file:
            elab_path = file
        elif output_identifier in file:
            output_path = file
    print('elab_path', elab_path)
    print('output_path', output_path)
    if elab_path == '' or output_path == '':
        print(f"Could not find elab or output csv file for {dir_path}")
        return None

def find_create_output_csv(root_path: str, dir_path: str):
    """
    Creates output.csv if it does not exist and there is an output folder. Should name it differently to avoid
    overwriting the original output.csv from Fragmenstein.

    INPUT:
        root_path: root path to parent directory
        dir_path: directory name
    """
    output_path = os.path.join(root_path, dir_path) + '/output'
    if not os.path.exists(output_path):
        print(f"{output_path} DOES NOT EXIST")
        return
    output_csv_path = os.path.join(root_path, dir_path, 'output.csv')
    if not os.path.exists(output_csv_path):
        print(f"{output_csv_path} DOES NOT EXIST. MAKING ONE NOW...")
        temp_output_csv_path = os.path.join(root_path, dir_path, 'temp_output.csv')
        # TODO: MAKE OUTPUT CSV SAME FORMAT AS FRAGMENSTEIN IN fragmenstein/laboratory/_base.py

        return



def get_merged_data(root_path: str, dir_path: str, elab_path: str, output_path: str,
                    rmsd_thresh: float) -> pd.DataFrame:
    """
    Gets dataframe of successful placements based on rmsd threshold.

    INPUT:
        dir_path: directory path of the base compound
        elab_path: path to the elab csv file
        output_path: path to the output csv file
    OUTPUT:
        acceptable_data: DataFrame of the acceptable placements and its metadata
        elab_path: path to the csv of elaborations
    """
    merged_data = []
    elabdf = pd.read_csv(elab_path)
    elabdf.drop_duplicates(inplace=True)  # always drop duplicates
    # Check if elabdf is already merged by if it contains 'outcome' column
    if 'outcome' not in elabdf.columns:  # Merge
        outputdf = pd.read_csv(output_path, index_col=0)
        elabdf = elabdf.merge(outputdf, on='smiles', how='left', suffixes=('_elab', '_output'))
    elabdf.to_csv(elab_path, index=False)
    assert 'comRMSD' in elabdf.columns, f"RMSD column does not exist! Please check {elab_path}"
    # Find delta delta G column
    # ddg_identifier = 'ΔΔG'
    ddg_identifier = unicodedata.normalize('NFKC', '∆∆G')
    ddg_column = None
    for col in elabdf.columns:
        col = unicodedata.normalize('NFKC', col)
        if ddg_identifier in col:
            ddg_column = col
    if ddg_column is None:
        print(f'∆∆G was not found! Please check {elab_path}...')
    acceptable_data = elabdf[(elabdf['comRMSD'] <= rmsd_thresh) & (elabdf[ddg_column] < 0)]
    return acceptable_data, elab_path


def move_folder(row, output_path, success_dir_path, df, remove: bool = True):
    """
    Moves the folder to the success dir. Also checks if the extra fragmenstein files are there, if so removes them.

    :param row:
    :param output_path:
    :param success_dir_path: path to success directory
    :param df:
    :param remove: boolean if the extra fragmenstein files should be removed

    :return:
    """
    # Check if 'moved_to_success_dir' column exists in df, if not, add it
    if 'moved_to_success_dir' not in df.columns:
        df['moved_to_success_dir'] = False
    # Define the source and destination paths
    source = os.path.join(output_path, row['name'].replace('_', '-'))
    destination = os.path.join(success_dir_path, row['name'].replace('_', '-'))
    # Check if the destination directory already exists
    if not os.path.exists(destination):
        # If it does not exist, attempt to move the folder
        try:
            shutil.move(source, destination)
            # Update the DataFrame if move is successful
            df.at[row.name, 'moved_to_success_dir'] = True
        except Exception as e:
            # If there's an error, print it, but continue with the next iteration
            print(f"Could not move {source} to {destination}: {e}")
    else:
        # If it exists, print a message and continue
        df.at[row.name, 'moved_to_success_dir'] = True
        print(f"Directory {destination} already exists, skipping.")
    # Check if the extra fragmenstein files are there, if so remove them
    if remove:
        suffixes_to_keep = ['.minimised.json', '.holo_minimised.pdb', '.minimised.mol']
        for root1, dirs, files in os.walk(success_dir_path):
            for directory2 in dirs:
                directory2_path = os.path.join(success_dir_path, directory2)
                if not os.path.exists(directory2_path):
                    print('THIS PATH DOES NOT EXIST', directory2_path)
                    continue
                for file in os.listdir(directory2_path):
                    if not any(file.endswith(suffix) for suffix in suffixes_to_keep):
                        file_path = os.path.join(success_dir_path, directory2, file)
                        try:
                            if os.path.isfile(file_path) or os.path.islink(file_path):
                                os.unlink(file_path)
                                print(f"Deleted file: {file_path}")
                            elif os.path.isdir(file_path):
                                shutil.rmtree(file_path)
                                print(f"Deleted directory: {file_path}")
                        except Exception as e:
                            print('Failed to delete %s. Reason: %s' % (file_path, e))


def create_success_directories(root: str, dir: str, acceptable_data: pd.DataFrame, elab_path: str, remove: bool) -> Tuple[
    Dict[str, int], int]:
    """
    Create success directories and return a mapping of paths to number of files. Will do an exhaustive check if success
    dirs already exist and will move the ones that have not been created yet.

    INPUT:
        root: root path to parent directory
        dir: directory name
        acceptable_data: dataframe of the metadata of the acceptable placements
        elab_path: path to the elaboration csv
    """
    success_dict = {}
    success_dir_path = os.path.join(root, dir, 'success')
    if not os.path.exists(success_dir_path):
        os.makedirs(success_dir_path)
    output_path = os.path.join(root, dir) + '/output'
    # move dirs in success.csv to success dir
    try:
        acceptable_data.apply(lambda row: move_folder(row, output_path=output_path, success_dir_path=success_dir_path,
                                                      df=acceptable_data, remove=remove), axis=1)
    except Exception as e:
        print(e)
    num_dirs = sum(os.path.isdir(os.path.join(success_dir_path, i)) for i in os.listdir(success_dir_path))
    success_dict[success_dir_path] = num_dirs
    acceptable_data.to_csv(os.path.join(success_dir_path, 'success.csv'))
    # Merge elabdf with acceptable data
    elabdf = pd.read_csv(elab_path)
    acceptable_data_subset = acceptable_data[['smiles', 'moved_to_success_dir']]
    elab_path2 = Path(elab_path)
    success_path_identifier = '_success_moved'
    # Check if 'moved_to_success_dir' already exists, if so, delete. Also delete suffix of
    # success_moved since adding it again will be redundant.
    if 'moved_to_success_dir' in elabdf.columns:
        del elabdf['moved_to_success_dir']
        # Remove the identifier from the stem
        new_stem = elab_path2.stem.replace(success_path_identifier, '')
        elab_path2 = elab_path2.with_name(new_stem + elab_path2.suffix)
    df = elabdf.merge(acceptable_data_subset, on='smiles', how='left', suffixes=('_elab', '_output'))
    df.drop_duplicates(inplace=True)
    success_dir_merge_path = elab_path2.with_stem(f"{elab_path2.stem}{success_path_identifier}")
    # Save merged data to new csv
    df.to_csv(success_dir_merge_path, index=False)
    # Save elab dataframe to same path
    elabdf.to_csv(elab_path, index=False)
    return success_dict, num_dirs


# ----------------------------#
def main():
    parser = config_parser()
    args = parser.parse_args()
    total_success_dirs = []
    num_success_dirs = 0
    num_success_elabs = 0
    for root, directories, files in os.walk(
            args.directory):  # walk through home directory, where the directories of the elab compounds are
        # Only interested in the top-level directories within args.directory
        if root != args.directory:
            # If the root is not the directory we started with, it's a subdirectory, so skip it
            continue
        for directory in directories:
            if directory == 'success' or directory == 'output' or 'xtra_results' in directory:
                continue
            # directory = name of the elab compound
            # Merge csvs of the original csv and output csv from Fragmenstein (if not already done)
            try:
                if args.rmsd is None:
                    rmsd_thresh = 2.0
                else:
                    rmsd_thresh = float(args.rmsd)
                # Create output.csv if it does not exist
                output_path = find_create_output_csv(root, directory, args.output, args.elab)
                elab_path = find_create_elab_csv(root, directory, args.output, args.elab)
                acceptable_data, elab_path = get_merged_data(root, directory, output_path, elab_path, rmsd_thresh)
                if acceptable_data is None:
                    print(f"{directory} DOES NOT CONTAIN THE ELAB AND OUTPUT CSV FILES")
                    continue  # go onto next directory of elab compounds
                # Create success dirs
                success_dict, num_dirs = create_success_directories(root, directory, acceptable_data, elab_path, args.remove)
                num_success_dirs += 1
                num_success_elabs += num_dirs
                # Store success dirs paths and number of successes within each
                total_success_dirs.append(success_dict)
                print(f"SUCCESS: {directory}")
            except Exception as e:
                print(f"Error in {directory}: {e}")
        break
    # Save results
    with open((os.path.join(args.directory, 'success_dirs.json')), 'w') as f:
        json.dump(total_success_dirs, f)


if __name__ == "__main__":
    main()
