#!/usr/bin/env python3

"""
Good lucking finding helpful analogues!

Author: Kate Fieseler
"""
import ast
import os
import sys
import argparse
import pandas as pd
from rdkit import Chem
import numpy as np
import glob2
import traceback
import datetime

sys.path.append('/Users/kate_fieseler/PycharmProjects/chemUtils')

from config import config
from wholeMoleculePipeline import searchReactantAnalogues, searchExactReactantAnalogues
from constants import REACTIONS_NAMES
from syndirella.utils import input


def find_rows_to_process(df, batch_num, batch_size):
    # Get rows to be processed based on batch_num, batch_size or row argument
    if args.row is not None:
        rows_to_process = [args.row]
    elif args.batch is not None:
        start_idx, end_idx = calculate_batch_indices(len(df), args.batch_size, args.batch)
        rows_to_process = range(start_idx, end_idx)
    else:
        rows_to_process = df.index
    return rows_to_process


def calculate_batch_indices(total_rows, batch_size, batch_num):
    """
    Calculate the start and end indices for a batch.
    """
    start_idx = (batch_num - 1) * batch_size
    end_idx = min(start_idx + batch_size, total_rows)
    return start_idx, end_idx


def searchAnalogues_2ndstep(df, results_dir, superstructure, rows_to_process=1, max_combinations=10_0000):
    """
    Given a dataframe of step 2 reactants with columns: smiles, num_steps, rxn_order_first_to_last, dir_name,
    reactant1_exact, reactant2_exact.
    Prepare input to searchReactantAnalogues()

    INPUT:
        :param df: step 2 reactants information.
        :param results_dir:
        :param superstructure:
        :param rows_to_process:
    OUTPUT:
        df: dataframe summarizing the number of analogues found for all the base compounds given in step2_df.
            Adding columns 'Analogs' and 'Done_Time' summarizing.
    """
    # Create output directory
    try:
        os.mkdir(results_dir)
    except IOError:
        pass
    if 'Done_Time' not in df.columns:
        df['Done_Time'] = np.nan
    if 'Analogs' not in df.columns:
        df['Analogs'] = np.nan
    # Convert NaN values to None
    df = df.applymap(lambda x: None if pd.isna(x) else x)

    for row_num in rows_to_process:
        row = df.iloc[row_num]
        if 'Done_Time' in df.columns:
            if row['Done_Time'] != None:
                continue
        if args.row is not None:
            if row_num != args.row:
                continue
        reaction_name = ast.literal_eval(row['rxn_order_first_to_last'])[0].replace(' ', '_')
        if reaction_name not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['smiles']}...\n")
            continue
        reactant1 = row['reactant1_exact']
        reactant2 = row['reactant2_exact']
        # TODO: Need to make sure that N-Bn deprotections are run since both reactants will not be found
        try:
            assert (reactant1 is None) != (reactant2 is None), ("Both reactants are provided in 2 step route, one must "
                                                                "be a product from step 1.")
        except AssertionError as e:
            print(e)
            continue
        ori_mol = Chem.MolFromSmiles((row['smiles']))
        which_elab = 0
        if reactant1 is not None:
            print('ELABORATING REACTANT 1')
            elab_reactant = reactant1
            which_elab = 1
        if reactant2 is not None:
            print('ELABORATING REACTANT 2')
            elab_reactant = reactant2
            which_elab = 2
        elab_reactant_mol = Chem.MolFromSmiles(elab_reactant)
        if elab_reactant_mol.GetNumAtoms() < 5:
            print(
                "Elaborated reactant has very little atoms, the full pipeline search will not be performed since that "
                "functionality is not implemented yet...\n")
            continue
        dir_names = [row['dir_name']]
        # Remove colon from dir names
        dir_names = [dir_name.replace(':', '_') for dir_name in dir_names]
        reaction_dir_names = []
        xtra_reaction_dir_names = []
        for i in range(len(dir_names)):  # put results in each directory
            dir_name = dir_names[i]
            reaction_dir_name = f"{results_dir}/{dir_name}/"
            os.makedirs(reaction_dir_name, exist_ok=True)
            reaction_dir_names.append(reaction_dir_name)
            xtra_reaction_dir_name = f"{results_dir}/{dir_name}/xtra_results_2/"
            os.makedirs(xtra_reaction_dir_name, exist_ok=True)
            xtra_reaction_dir_names.append(xtra_reaction_dir_name)
            output_name = f"{dir_name}_2nd_step"
            # find product_df, in first level of reaction_dir_name, assert it is NOT output.csv
            csv_files = glob2.glob(f'{reaction_dir_name}*.csv')
            product_path = None
            for file in csv_files:
                if dir_name in file and 'output' not in file and '2nd' not in file and '2_of_2' not in file:
                    product_path = file
            if product_path is None:
                print(f'NO PRODUCT CSV FOUND IN {dir_name}')
                continue
            product_df = pd.read_csv(product_path)
            assert which_elab != 0, "The elaborated reactant has not been identified! Stopping..."
            output_name = f"{dir_name}_2nd_step"  # TODO: change this to be more descriptive
        try:
            results = searchExactReactantAnalogues(elab_reactant_mol, product_df, which_elab, original_mol=ori_mol,
                                                   resultsDirs=xtra_reaction_dir_names,
                                                   reaction_name=reaction_name,
                                                   similarityThr=config.SIMILARITY_SEARCH_THR,
                                                   structuralThr=config.STRUCTURAL_SCORE_THR, output_name=output_name,
                                                   max_combinations=max_combinations)
            # Set the Done_Time and Analogs columns for the current row
            df.at[row_num, 'Done_Time'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            df.at[row_num, 'Analogs'] = len(results) if results is not None else 0
        except Exception as e:
            tb = traceback.format_exc()
            print(f"Error processing row {i}: {tb}")
            continue
        if results is None:
            print("No results found for this molecule.\n")
            continue
        for i in range(len(reaction_dir_names)):
            num_analogs = len(results)
            print('SAVING RESULTS AT:', reaction_dir_names[i])
            results.to_csv(os.path.join(reaction_dir_names[i], f"{output_name}_{num_analogs}_2_of_2.csv"), index=False)
        print(results)
    if "OUTPUT" in args.input_csv:
        df.to_csv(args.input_csv, index=True)
    else:
        df.to_csv(f'{args.input_csv.split(".")[0]}_OUTPUT.csv', index=False)


def searchAnalogues(df, results_dir, superstructure, step_num, rows_to_process=1):
    '''
    :param df: dataframe with columns SMILES, Reaction_name, Reactants
    :param results_dir: directory to save results
    :param superstructure: if True, perform superstructure search
    :return:
    '''
    # Create output directory
    try:
        os.mkdir(results_dir)
    except IOError:
        pass
    if 'Done_Time' not in df.columns:
        df['Done_Time'] = np.nan
    if 'Analogs' not in df.columns:
        df['Analogs'] = np.nan
    # Convert NaN values to None
    df = df.applymap(lambda x: None if pd.isna(x) else x)
    assert df['smiles'] is not None, "No SMILES provided in input csv."

    for row_num in rows_to_process:
        start_time = datetime.datetime.now()  # Start timing
        row = df.iloc[row_num]
        if 'Done_Time' in df.columns:
            if row['Done_Time'] != None:
                continue
        if args.row is not None:
            if row_num != args.row:
                continue
        step = row['num_steps']
        """
        This part makes directories for each route and csv and runs the whole molecule pipeline
        
        INPUT TO WHOLE MOLECULE PIPELINE:
            ori_mol: RDKit mol object of the product
            reactant1_mol: RDKit mol object of the first reactant
            reactant2_mol: RDKit mol object of the second reactant
            ori_reaction: Reaction name to produce product
            resultsDirs: Path to save xtra results
            output_name: Name of the output csv (not full path)
        """
        ori_mol = Chem.MolFromSmiles((row['smiles']))
        if row['num_steps'] == 1:
            reaction_name = ast.literal_eval(row['rxn_order_first_to_last'])[0].replace(' ', '_')
        if reaction_name not in REACTIONS_NAMES:
            print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
                  f"Please provide the SMARTS.\n"
                  f"Skipping {row['smiles']}...\n")
            continue
        else:
            reactants = ast.literal_eval(row['reactants'])[0]  # CHANGE FOR OTHER STEPS
            print(reactants[0])
            reactant1_mol = Chem.MolFromSmiles(reactants[0])
            reactant2_mol = Chem.MolFromSmiles(reactants[1])
            if reactant1_mol.GetNumAtoms() < 5 or reactant2_mol.GetNumAtoms() < 5:
                print("One reactant has very little atoms, the full pipeline search will not be performed since that "
                      "functionality is not implemented yet...\n")
                continue
            dir_names = [row['dir_name']]
            # Remove colon from dir names
            dir_names = [dir_name.replace(':', '_') for dir_name in dir_names]
            reaction_dir_names = []
            xtra_reaction_dir_names = []

            # Get constant fragmenstein placement info to append in dictionary
            fragmen_info = {'hit_names': row.hit_names, 'ref_pdb': row.ref_pdb}

            for i in range(len(dir_names)):
                dir_name = dir_names[i]
                reaction_dir_name = f"{results_dir}/{dir_name}/"
                os.makedirs(reaction_dir_name, exist_ok=True)
                reaction_dir_names.append(reaction_dir_name)
                xtra_reaction_dir_name = f"{results_dir}/{dir_name}/xtra_results_{step}/"
                os.makedirs(xtra_reaction_dir_name, exist_ok=True)
                xtra_reaction_dir_names.append(xtra_reaction_dir_name)
                output_name = f"{dir_name}_{step}_of_1_step"
            try:
                results = searchReactantAnalogues(ori_mol, reactant1_mol, reactant2_mol, step_num=step,
                                                  ori_reaction=reaction_name,
                                                  resultsDirs=xtra_reaction_dir_names, output_name=output_name,
                                                  struct_score=False, fragmen_info=fragmen_info)
                # Set the Done_Time and Analogs columns for the current row
                df.at[row_num, 'Done_Time'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                df.at[row_num, 'Analogs'] = len(results) if results is not None else 0
            except Exception as e:
                tb = traceback.format_exc()
                print(f"Error processing row {i}: {tb}")
                num_analogs = 0
                continue
            if results is None:
                print("No results found for this molecule.\n")
                continue
            for i in range(len(reaction_dir_names)):
                num_analogs = len(results)
                results.to_csv(os.path.join(reaction_dir_names[i], f"{output_name}_{num_analogs}_1_of_1.csv"),
                               index=False)
            print(results)
            end_time = datetime.datetime.now()  # End timing
            duration = end_time - start_time  # Calculate duration
            print(f"row {i} took {duration} to execute")
    if "OUTPUT" in args.input_csv:
        df.to_csv(args.input_csv, index=True)
    else:
        df.to_csv(f'{args.input_csv.split(".")[0]}_OUTPUT.csv', index=False)

    return os.path.join(reaction_dir_names[i], f"{output_name}_{num_analogs}.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieve synthesizable analogues from csv file of product SMILES, reactants, and reaction name.")
    parser.add_argument('-i', '--input_csv', type=str,
                        help=('Path to the input CSV file. The expected CSV structure is:\n'
                              'SMILES (str) - SMILES of product\n'
                              'dir_name (list) - Name of the directories to save results. Usually target id.\n'
                              'num_steps (int) - Number of steps in the route\n'
                              'rxn_order_first_to_last (list(str)) - Reaction name to produce product\n'
                              'reactants (list(tuple)) - Reactants listed in tuples\n'
                              '...\n'))
    parser.add_argument('-r', "--results_dir", help="Directory for the results", required=True)
    parser.add_argument('-u', "--superstructure", help='if performing a superstructure search', action="store_true")
    parser.add_argument('-b', '--row', help='specify row number to search. 0 is the first row below the header.',
                        type=int)
    parser.add_argument('--batch', help='specify batch number to search. 1 is the first batch.', type=int)
    parser.add_argument('--batch_size', help='specify batch size. Default is 10.', type=int, default=10)
    parser.add_argument('--go', help='Specify if you want to keep going after starting a batch.', action="store_true")
    parser.add_argument('-n', '--num_steps', help='Specify the number of steps in the route.', type=int, default=1)
    parser.add_argument('--test', help='Specify if you want to test the script.', action="store_true")
    parser.add_argument('--step1',
                        help=".csv of the first step of multi step routes. Only needed for more than 1 step route.",
                        type=str)
    parser.add_argument('--step2',
                        help=".csv of the second step of multi step route. Only needed for more than 1 step route.",
                        type=str)
    parser.add_argument('--max_combinations', help='Specify the maximum number of combinations to search for.',
                        type=int, default=10000)

    args = parser.parse_args()

    # Just run this for elaborating step 2
    if args.test and args.step2:
        step2_df = pd.read_csv(args.step2, sep=',', header=0)
        rows_to_process = find_rows_to_process(step2_df, args.batch, args.batch_size)
        searchAnalogues_2ndstep(step2_df, args.results_dir, args.superstructure, rows_to_process=rows_to_process,
                                max_combinations=args.max_combinations)
        exit()

    # if args.test and args.num_steps==2:
    #     # Load the smiles into dataframe
    #     df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
    #     # Get rows to be processed based on batch_num, batch_size or row argument
    #     rows_to_process = find_rows_to_process(df, args.batch, args.batch_size)
    #     searchAnalogues(df, args.results_dir, args.superstructure, step_num=2, rows_to_process=rows_to_process)

    if args.test and args.num_steps == 2:
        df = pd.read_csv(args.input_csv, sep=',', header=0)
        if args.step1 is not None:
            step1_df = pd.read_csv(args.step1, sep=',', header=0, index_col=0)
            if args.step2 is not None:
                step2_df = pd.read_csv(args.step2, sep=',', header=0, index_col=0)
            else:
                step2_df = input.editstep2(df, step1_df)
                step2_df.to_csv(f'{args.input_csv.split(".")[0]}_2nd_step_TEST.csv', index=True)
        else:
            df = pd.read_csv(args.input_csv, sep=',', header=0)
            step1_df = input.editstep1(df)
            step2_df = input.editstep2(df, step1_df)
            step1_df.to_csv(f'{args.input_csv.split(".")[0]}_1st_step_TEST.csv', index=False)
            step2_df.to_csv(f'{args.input_csv.split(".")[0]}_2nd_step_TEST.csv', index=False)
        # df_edit = editstep1forstep2(df)
        # searchAnalogues(df_edit, args.results_dir, args.superstructure, step_num=2)

    # Load the smiles into dataframe

    df = pd.read_csv(args.input_csv, sep=',', header=0)
    # Get rows to be processed based on batch_num, batch_size or row argument
    rows_to_process = find_rows_to_process(df, args.batch, args.batch_size)
    searchAnalogues(df, args.results_dir, args.superstructure, step_num=1, rows_to_process=rows_to_process)

    if args.num_steps == 2:
        df = pd.read_csv(args.input_csv, sep=',', header=0, index_col=0)
        df = input.editstep1(df)
        df.to_csv(f'{args.input_csv}_1st_step', index=True)
        rows_to_process = find_rows_to_process(df, args.batch, args.batch_size)
        step1_output = searchAnalogues(df, args.results_dir, args.superstructure, step_num=1)
        df = pd.read_csv(step1_output, sep=',', header=0, index_col=0)
        df_edit = editstep1forstep2(df)
        searchAnalogues(df_edit, args.results_dir, args.superstructure, step_num=2)
