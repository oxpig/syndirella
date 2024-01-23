#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:12:49 2021

@author: ruben
"""
import itertools
import os
from collections import OrderedDict
from typing import Tuple, Iterable
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

import numpy as np
import pandas as pd
from chemUtils.geometry.molDistances import findCloseAtomIdxsFromRefAtomIdxs, find_exact_match_atom_indices, find_substruct_atom_match
from chemUtils.substructure import splitExtendedMolBySubstruct
from chemUtils.substructure import find_attachmentIdxs_fromMolAndSubstruct
from rdkit.Chem import rdFMCS
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.PandasTools import WriteSDF
from rdkit.rdBase import BlockLogs

from config import config
from constants import fromReactionFullNameToReactantName, REACTIONS_NAMES
from databaseSearch.databaseSearch import perform_database_search
from embed import embedMolUsingRefs
from engine import myMap
from reactionSmarts import (checkReactionSmartInAttachment, PATTERN_REACTANT1, PATTERN_REACTANT2, \
    runReactionInAttachment, checkReactionSmartInAllAtoms, reaction_smarts,
                            checkReactionSmartInAllAtomsAndReactants, checkReactionSmartAroundAttachment, checkSpecificReactionSmartInReactant,
                            REACTANT1_DICT, REACTANT2_DICT)
from scorer import scorePairOfMols

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdfiltercatalog import FilterCatalog, FilterCatalogParams

import pickle
from itertools import chain

from PIL import Image, ImageDraw, ImageFont
import ast

def add_fragmenstein_name(row, output_name, column_names):
    """
    Adds the fragmenstein name to the dataframe

    :param row: A pandas Series representing a row of the DataFrame
    :param output_name: A string to be used in the naming
    :param column_names: A list of column names in the DataFrame
    :return: Updated row with the name added
    """
    if 'reactant1_structuralScore' in column_names and 'reactant2_structuralScore' in column_names:
        if row['reactant1_structuralScore'] == 1 and row['reactant2_structuralScore'] == 1:
            row['name'] = f"{output_name}-base"
    else:
        structuralScore_cols = [col for col in column_names if 'structuralScore' in col]
        if structuralScore_cols and 'num_atom_difference' in column_names:
            structuralScore_col = structuralScore_cols[0]  # Assuming you want to use the first match
            if row[structuralScore_col] == 1 and row['num_atom_difference'] == 0:
                row['name'] = f"{output_name}-base"
            else:
                row['name'] = f"{output_name}-{row['index']}"
    return row

def atom_num_difference(base, e_mol):
    mcs = rdFMCS.FindMCS([base, e_mol])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    mcs_atoms = mcs_mol.GetNumAtoms()
    e_mol_atoms = e_mol.GetNumAtoms()
    difference = e_mol_atoms - mcs_atoms
    return difference

def pains_filter(df):
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    def is_pain(smile):
        if not isinstance(smile, str):
            return None
        mol = Chem.MolFromSmiles(smile)
        if mol and catalog.HasMatch(mol):
            return True
        return False

    df['is_pains'] = df['smiles'].apply(is_pain)

    return df

def save_molecule_with_indices(mol, filename, text):
    # Generate a new molecule that includes atom indices
    mol_with_indices = Chem.Mol(mol)
    for atom in mol_with_indices.GetAtoms():
        atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))

    # Convert the molecule with indices to an image
    img = Draw.MolToImage(mol_with_indices, size=(500, 500))

    # Convert RDKit Image to PIL Image
    pil_img = Image.new("RGB", img.size)
    pil_img.paste(img)

    # If there's text to add
    if text:
        draw = ImageDraw.Draw(pil_img)
        draw.text((10, 10), text, (255, 255, 255))  # Coordinates, text, color, and font

    # Save the image to the file
    pil_img.save(filename)

def getFunction_checkStructuralMatch(originalMol):
    originalMol = Chem.AddHs(originalMol, addCoords=True)

    def analyzeOneSmi(simMol):
        # if not simMol:
        #     return None
        logsManager = BlockLogs()
        # Changing mol to explicit Hs which should help with embedding
        simMol = Chem.AddHs(simMol, addCoords=True)
        simEmbedded = embedMolUsingRefs(simMol, originalMol)  # TODO: fix to not return None all the time
        AllChem.EmbedMolecule(originalMol) # Let_s try to embed the original mol, it worked
        # calculate the mw of original and new mol
        mw_original = Descriptors.MolWt(originalMol)
        mw_sim = Descriptors.MolWt(simMol)
        mw_diff = abs(mw_original - mw_sim)
        del logsManager
        if Chem.MolToSmiles(originalMol) == Chem.MolToSmiles(simMol): # keeps original reactant
            return 1, originalMol, mw_diff, mw_sim
        if simEmbedded is None:
            return None, None, None, None
        else:
            return scorePairOfMols(originalMol, simEmbedded)[0], simEmbedded, mw_diff, mw_sim

    return analyzeOneSmi


def embedAndScoreSimilarSmiles(refMol: Chem.Mol, smilesList) -> Tuple[Iterable[float], Iterable[Chem.Mol]]:
    """

    :param refMol: The reference mol to embed the smiles
    :param smilesList: A list of smiles
    :return: a list of tuples (structScore, smiles, similarity)
    """

    similar_mols = myMap(Chem.MolFromSmiles, smilesList)
    results = myMap(getFunction_checkStructuralMatch(refMol), similar_mols)

    structuralScore, embedded_similar_mols, mw_diff, mw = zip(*results) if results else ([], [], [], [])
    return structuralScore, embedded_similar_mols, mw_diff, mw


def matchMolToReactan(df, reaction_name, reactant_smarts, reactant_Num):
    """
    Given a dataframe with a column "mol" and a column "attachmentIdx", this function checks if the molecule is a
    reactant of specific reaction provided with reaction_name.

    :param df:
    :param reactantNum:
    :return:
    """
    # if reactantNum == 1:
    #     valid_patterns = pattern_reactant1
    # elif reactantNum == 2:
    #     valid_patterns = pattern_reactant2  # Need to have two reactants, flip the order to find both combinations
    # else:
    #     raise ValueError("Only 1 and 2 valid for reactant_num")
    # matching_reactions = myMap(lambda mol_attachmentIdx:
    #                            checkReactionSmartAroundAttachment(*mol_attachmentIdx,
    #                                                           valid_patterns=valid_patterns, reactantNum=reactantNum),
    #                            zip(df["mol"], df["attachmentIdx"]))
    matching_reactions = myMap(lambda smiles:
                               checkSpecificReactionSmartInReactant(smiles, reaction_name, reactant_smarts),
                               df["smi"])

    selected_reactions = {name: [] for name in REACTIONS_NAMES}
    # can return false here
    for i, reaction_name_flag in enumerate(matching_reactions):
        for name, flag in reaction_name_flag:
            selected_reactions[name].append(flag)

    vals = selected_reactions[reaction_name]
    assert len(vals) == len(df), ("Error, the number of values in the selected reactions is not the same as the number "
                                  "of molecules")
    df[fromReactionFullNameToReactantName(name, reactant_Num)] = vals  # making new column in dataframe

    return df


def findReactionAtoms(ori_mol, reactant, reaction_name, results_dir, react_num, attachmentIdxs):
    """
    Given a molecule and a reactant, this function finds the atoms that are involved in the reaction for the reactant.
    :param ori_mol:
    :param reactant:
    :param reaction_name:
    :return:
    """
    matched_reaction_atoms = checkReactionSmartInAllAtoms(reactant, reaction_name, )
    if matched_reaction_atoms is None:
        print('No reaction atoms found for reactant', Chem.MolToSmiles(reactant), 'in reaction', reaction_name)
        return None
    matched_reaction_atoms = tuple(x for sub in matched_reaction_atoms for x in sub)  # This flattens the tuple
    save_molecule_with_indices(reactant, os.path.join(results_dir,
                                                      f"{reaction_name}_reactant{react_num}.png"), matched_reaction_atoms)

    # Check reaction atoms are in attachment idx
    if not set(attachmentIdxs).issubset(set(matched_reaction_atoms)): # attachmentIdxs is always smaller than matched_reaction_atoms (should be?)
        print('Reaction atoms are not in attachment idxs')
        return None
    else:
        # Intersection of the attachmentIdxs and the matched_reaction_atoms, don't neccesarily want this
        if reaction_name == 'Amidation' and len(matched_reaction_atoms) == 3:
            return matched_reaction_atoms
        matched_reaction_atoms = list(set(attachmentIdxs).intersection(set(matched_reaction_atoms)))
    return matched_reaction_atoms

def findReactionAtoms_bothReactants(originalMol, reactant1, reactant2, ori_reaction, resultsDirs, reactant1_attachmentIdxs, reactant2_attachmentIdxs):
    """
    Given a product, two reactants and a labeled reaction between them, need to find the atoms involved in the reaction on the reactants.

    :param originalMol:
    :param reactant1:
    :param reactant2:
    :param ori_reaction:
    :param resultsDirs: list of results dirs
    :param reactant1_attachmentIdxs:
    :param reactant2_attachmentIdxs:
    :return:
    """

    # matched_reaction_atoms is a dictionary of 'reactant1' and 'reactant2' with values as atom indices involved in reaction
    matched_reaction_dict = checkReactionSmartInAllAtomsAndReactants(reactant1, reactant1_attachmentIdxs, reactant2, reactant2_attachmentIdxs, reaction_name=ori_reaction)
    matched_reaction_atoms1 = matched_reaction_dict['reactant1'][2]
    matched_reaction_atoms2 = matched_reaction_dict['reactant2'][2]
    # Need to save with correct image label corresponding to SMARTS order
    reactant1 = matched_reaction_dict['reactant1'][1]
    reactant2 = matched_reaction_dict['reactant2'][1]
    for i in range(len(resultsDirs)):
        save_molecule_with_indices(reactant1, os.path.join(resultsDirs[i],
                                                          f"{ori_reaction}_reactant{1}.png"), matched_reaction_atoms1)
        save_molecule_with_indices(reactant2, os.path.join(resultsDirs[i],
                                                           f"{ori_reaction}_reactant{2}.png"), matched_reaction_atoms2)

    matched_reaction_info = ((matched_reaction_dict['reactant1'][0], matched_reaction_dict['reactant1'][1], matched_reaction_dict['reactant1'][2]), (matched_reaction_dict['reactant2'][0], matched_reaction_dict['reactant2'][1], matched_reaction_dict['reactant2'][2]))
    return matched_reaction_info


def createAnaloguesDf(smisDict, structuralScores, mols, mw_diff, mw, similarity, reactant):
    if similarity:
        similarities, metadatas = zip(*smisDict.values())
        result = OrderedDict(structuralScore=structuralScores, mw_diff=mw_diff, mw=mw, similarities=similarities, smi=list(smisDict.keys()),
                             metadata=metadatas, mol=mols)
    else:
        metadatas = list(smisDict.values())
        result = OrderedDict(structuralScore=structuralScores, mw_diff=mw_diff, mw=mw, smi=list(smisDict.keys()),
                             metadata=metadatas, mol=mols)
    result = pd.DataFrame(result)
    result.dropna(axis=0, inplace=True)
    return result


def processSubMolReactant(reactant, reaction_name, reactant_smarts, reaction_atoms, similarityThr, structuralThr, reactant_Num, resultsDirs, substructure=True):
    """
    Given a reactant and the atoms involved in the reaction, this function searches for analogues of the reactant
    with Postera's superstructure search.

    :param reactant: mol object of reactant
    :param reaction_atoms: atoms involved in the reaction
    :param similarityThr:
    :param structuralThr:
    :param resultsDirs: list of results dirs
    :param substructure: True if searching for superstructure
    :return: DataFrame with analogues of the reactant.
    """
    print('This is the reactant:', Chem.MolToSmiles(reactant))
    print('These are the atoms involved in the reaction:', reaction_atoms)

    # we know which reactant it is
    smisDict = perform_database_search(reactant, substructure=substructure, thr=similarityThr) # substructure is the search for superstructure

    # Save smisDict to a file, to show structures easily
    for i in range(len(resultsDirs)):
        with open(os.path.join(resultsDirs[i], 'smisDict.pickle'), 'wb') as f:
            pickle.dump(smisDict, f)

    structuralScores, similarMols, mw_diff, mw = embedAndScoreSimilarSmiles(reactant, smisDict.keys())

    similarity = True
    if substructure is True: similarity = False

    df = createAnaloguesDf(smisDict, structuralScores, similarMols, mw_diff, mw, similarity, reactant)
    for i in range(len(resultsDirs)):
        df.to_csv(os.path.join(resultsDirs[i], f'analogues_preprocess_{Chem.MolToSmiles(reactant)}.csv'))

    # Search through elaborated compounds for the reaction atoms as substructure search
    #equiv_attachmentIdxs = myMap(lambda mol: find_substruct_atom_match(mol, reactant, reaction_atoms), df["mol"])

    # Have to create conformers for the elaborated mols and reactants
    # Do you need a 3D conformer? I'll make a 3D one. I think the 3D conformer was already calcuated when they were embedded to compare similarity


    # equiv_attachmentIdxs2 =  myMap(lambda mol: find_exact_match_atom_indices(mol, reactant, reaction_atoms, test=False), df['mol'])
    #
    # # The equiv_attachementIdx2 is
    #
    # df.loc[:, "attachmentIdx"] = equiv_attachmentIdxs2
    # df = df[df["attachmentIdx"].notna()]
    print('This is the length of the dataframe after substructure search:', len(df))
    if structuralThr is not None:
        df = df[df["structuralScore"] > structuralThr]
    df = matchMolToReactan(df, reaction_name, reactant_smarts, reactant_Num)  # should be specific reactant
    df = df.reset_index(drop=True)
    for i in range(len(resultsDirs)):
        df.to_csv(os.path.join(resultsDirs[i], f'analogues_substruct_{Chem.MolToSmiles(reactant)}_{len(df)}.csv'))

    return df

import pandas as pd
from itertools import product

def apply_reaction(row, original_mol: Chem.Mol, reaction_name) -> list():
    """
    Applies a reaction to a row of a dataframe. Need to make sure that the reaction SMARTS is flipped so both sides of
    the reaction are tried out. But need to make sure that the reactants are in the right order like how they were originally discovered.
    Since functional groups decorating a ring do not want to be reacted together even if they match the smarts pattern.

    :param row:
    :param reaction_name:
    :return: list of product SMILES
    """
    # Reactant1 is just the first reactant in the SMARTS, and reactant2 is the second reactant in the SMARTS
    smarts_pattern = reaction_smarts[reaction_name]
    reaction = AllChem.ReactionFromSmarts(smarts_pattern)

    reactant1 = Chem.MolFromSmiles(row['smi_reactant1'])
    reactant2 = Chem.MolFromSmiles(row['smi_reactant2'])

    products = reaction.RunReactants((reactant1, reactant2))

    # Check the products can be sanitized
    def can_be_sanitized(mol):
        try:
            Chem.SanitizeMol(mol)
            return True
        except:
            return False

    sani_products = [product[0] for product in products if can_be_sanitized(product[0])]

    # Convert the products to SMILES strings
    product_smiles = [Chem.MolToSmiles(product) for product in sani_products]

    if len(product_smiles) == 0:
        print(f'NO PRODUCTS FOUND FOR {row["smi_reactant1"]} and {row["smi_reactant2"]}')

    # Check for duplicate smiles, remove them
    product_smiles = list(set(product_smiles))
    atom_difference = [atom_num_difference(original_mol, Chem.MolFromSmiles(product)) for product in product_smiles]

    row['smiles'] = product_smiles
    row['num_atom_difference'] = atom_difference

    return row

def find_stereoisomers(smiles) -> list():
    # This function should return a list of stereoisomers for the given SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    # Generate all stereoisomers
    opts = StereoEnumerationOptions(unique=True)
    isomers = list(EnumerateStereoisomers(mol, options=opts))
    isomer_list = [Chem.MolToSmiles(isomer, isomericSmiles=True) for isomer in isomers]
    return isomer_list

def enumerate_stereoisomers(df):
    """
    Given a DataFrame with a column of smiles strings, this function enumerates all the stereoisomers of the smiles.
    it also needs a names column to add conformer letters to

    :param df:
    :return:
    """
    assert 'smiles' in df.columns, 'No smiles column in DataFrame'
    assert 'name' in df.columns, 'No name column in DataFrame'

    new_rows = []

    for index, row in df.iterrows():
        stereoisomers = find_stereoisomers(row['smiles'])

        for i, iso in enumerate(stereoisomers):
            new_row = row.copy()
            new_row['smiles'] = iso
            new_row['name'] = f"{row['name']}_{chr(65 + i)}"  # Appending A, B, C, etc., to the name
            new_row['conformer'] = chr(65 + i)
            new_rows.append(new_row)

    new_df = pd.DataFrame(new_rows)
    exploded_df = pd.concat([df, new_df], ignore_index=True)
    # remove NaNs
    exploded_df = exploded_df.dropna()
    return exploded_df

def findReactionCombinations(original_mol, ori_reaction, analogs_reactant1, analogs_reactant2, output_name: str, max_combinations: int = None):
    """
    Given a reaction and the analogues of the reactants, this function finds the reaction combinations that are possible.

    :param ori_reaction:
    :param analogs_reactant1:
    :param analogs_reactant2:
    :param resultsDir:
    :return: DataFrame with reactant combinations
    """
    i = ori_reaction

    column_name_1 = fromReactionFullNameToReactantName(i, 1)
    column_name_2 = fromReactionFullNameToReactantName(i, 2)

    # Function to filter DataFrame based on column existence and value
    def filter_df(df, column_name):
        if column_name in df.columns:
            return df[df[column_name] == True]
        else:
            return pd.DataFrame(columns=df.columns)

    # Filter both DataFrames based on the react1_ReactionName and react2_ReactionName columns
    filtered_df1 = filter_df(analogs_reactant1, column_name_1)
    filtered_df2 = filter_df(analogs_reactant2, column_name_2)

    # Switching reactant numbers
    if filtered_df2.empty and filtered_df1.empty:
        filtered_df2 = filter_df(analogs_reactant1, column_name_2)
        filtered_df1 = filter_df(analogs_reactant2, column_name_1)

    # Get the Cartesian product of the filtered DataFrame based on index
    reactant_combinations = pd.MultiIndex.from_product([filtered_df1.index, filtered_df2.index], names=['reactant1', 'reactant2']).to_frame(index=False)

    # Filter the columns to keep
    columns_to_keep = ['structuralScore', 'mw', 'smi', 'num_atom_difference']
    metadata_columns_df1 = filtered_df1.columns[filtered_df1.columns.str.contains('metadata')]
    metadata_columns_df2 = filtered_df2.columns[filtered_df2.columns.str.contains('metadata')]
    # Check if columns exist in filtered_df1 and filter the existing ones
    filtered_df1 = filtered_df1.loc[:,
                   filtered_df1.columns.isin(columns_to_keep) | filtered_df1.columns.isin(metadata_columns_df1)]
    # Check if columns exist in filtered_df2 and filter the existing ones
    filtered_df2 = filtered_df2.loc[:,
                   filtered_df2.columns.isin(columns_to_keep) | filtered_df2.columns.isin(metadata_columns_df2)]

    r1_percent_decrease = round((((len(analogs_reactant1)-len(filtered_df1))/(len(analogs_reactant1)))*100),2)
    r2_percent_decrease = round((((len(analogs_reactant2)-len(filtered_df2))/(len(analogs_reactant2)))*100),2)
    print(f'NUMBER OF REACTABLE UNIQUE ELABORATIONS OF REACTANT 1: {len(filtered_df1)}. {r1_percent_decrease}% decrease.')
    print(f'NUMBER OF REACTABLE UNIQUE ELABORATIONS OF REACTANT 2: {len(filtered_df2)}. {r2_percent_decrease}% decrease.')

    filtered_df1 = filtered_df1.add_suffix('_reactant1')
    filtered_df2 = filtered_df2.add_suffix('_reactant2')

    # Merge the information based on the indices
    result_df = (
        reactant_combinations
        .merge(filtered_df1, left_on='reactant1', right_index=True)
        .merge(filtered_df2, left_on='reactant2', right_index=True)
    )
    result_df = result_df.drop(columns=['reactant1', 'reactant2'])

    if max_combinations is not None and len(result_df) > max_combinations:
        print(f'Limiting number of reactant combinations to {max_combinations}. Randomly sampling.')
        # Separate the first 5 rows
        first_five_rows = result_df.head(5)
        # If the DataFrame has more than 5 rows, sample from the rest of the DataFrame
        if len(result_df) > 5:
            remaining_rows = result_df.iloc[5:]
            sampled_rows = remaining_rows.sample(n=max_combinations - 5, random_state=1)
            result_df = pd.concat([first_five_rows, sampled_rows])
        else:
            result_df = first_five_rows

    # Add a column of the reacted together SMILES
    result_df = result_df.apply(apply_reaction, original_mol=original_mol, reaction_name=i, axis=1)
    # Use the explode method to transform lists in the 'product_smi' column into separate rows
    result_df = result_df.explode(['smiles','num_atom_difference']).reset_index(drop=True)

    # Remove empty rows with no products
    result_df2 = result_df.copy()
    result_df2 = result_df2[result_df2['smiles'].notna()]
    print(f'{len(result_df) - len(result_df2)} empty product routes removed.')

    # Remove repeated rows
    result_df2 = result_df2.drop_duplicates(subset=['smi_reactant1', 'smi_reactant2', 'smiles'])
    print(f'{len(result_df)-len(result_df2)} product duplicates removed.')

    # How many unique products are there?
    print(f'{result_df2["smiles"].nunique()} unique products with unique routes found.')

    # SORT BY NUMBER OF ATOM DIFFERENCE
    result_df2.sort_values(by=['num_atom_difference'], inplace=True, ignore_index=True)

    # ADDING PAINS FILTER
    # pains_filter(result_df2)

    # ADD NAME FOR FRAGMENSTEIN PLACEMENT
    column_names = result_df2.columns.tolist()
    result_df2['index'] = result_df2.index
    result_df2 = result_df2.apply(add_fragmenstein_name, output_name=output_name, column_names=column_names, axis=1)
    # Check if multiple products have the same name
    # If so, add a number to the end of the name
    result_df2['count'] = result_df2.groupby('name').cumcount() + 1
    result_df2['name'] = result_df2.apply(lambda x: f"{x['name']}_{x['count']}" if x['count'] > 1 else x['name'],
                                          axis=1)
    result_df2.drop(columns='count', inplace=True)

    # Only take the top 10K products to enumerate stereoisomers
    if len(result_df2) > 10_000:
        print('Limiting number of products to 10,000.')
        result_df2 = result_df2.head(10_000)

    print('Enumerating stereoisomers...')
    result_df3 = enumerate_stereoisomers(result_df2)
    all_products = result_df3['smiles']

    # Save the results to a CSV file
    return result_df3, all_products

def reactBackTogether(original_mol: Chem.Mol, analogs_reactant1: pd.DataFrame, analogs_reactant2: pd.DataFrame, ori_reaction, output_name: str, max_combinations: int = None):
    if ori_reaction is None:
        print('Must provide reaction to react back together.')
        # TODO: Could provide a way to react back together with all other reactions in the database.
        return None
    else:
        df_combined_results, all_products = findReactionCombinations(original_mol, ori_reaction, analogs_reactant1, analogs_reactant2, output_name=output_name, max_combinations=max_combinations)

    return df_combined_results

def productsToReact(df, reactant_num, reactant_smarts, reaction_name, output_name='test'):
    """
    Given a dataframe of products, convert it to a dataframe of reactants for the second step.

    INPUT:
        df:
        reactant_num: int of reactant in SMARTS

    OUTPUT:
        product_react_df:
    """
    assert reactant_num is not None, 'Must provide reactant number to convert products to reactants.'
    # Function to concatenate JSON data
    def concat_json(row):
        reactant1: list = ast.literal_eval(row['metadata_reactant1'])
        reactant2: list = ast.literal_eval(row['metadata_reactant2'])
        return pd.Series({'metadata': [reactant1, reactant2]})
    def check_reactant_smarts(row):
        matched_reaction = checkSpecificReactionSmartInReactant(row['smi'], reaction_name, reactant_smarts)[0]
        if matched_reaction[1] is None or matched_reaction[1] == False:
            row[f'react{reactant_num}_{reaction_name}'] = False
        else:
            row[f'react{reactant_num}_{reaction_name}'] = True
        return row
    # Creating a new DataFrame with the required columns
    product_react_df = pd.DataFrame({
        'smi': df['smiles'],
        'num_atom_difference': df['num_atom_difference'],
        'metadata_reactant1': df['metadata_reactant1'],
        'metadata_reactant2': df['metadata_reactant2']
    })
    # Concatenate metadata
    product_react_df = pd.concat([product_react_df, df.apply(concat_json, axis=1)], axis=1)
    # Check SMARTS of reactant is in product
    product_react_df = product_react_df.apply(check_reactant_smarts, axis=1)
    # change index
    product_react_df.reset_index(drop=True, inplace=True)
    return product_react_df

def searchExactReactantAnalogues(elab_reactant, product_df, which_elab, original_mol, resultsDirs,
                                 reaction_name=None, similarityThr=config.SIMILARITY_SEARCH_THR,
                                 structuralThr=config.STRUCTURAL_SCORE_THR, output_name='test',
                                 max_combinations=10_000):
    """
    Given two exact reactants for a given reaction, elaborate one of them then find the product of those elaborations
    with the other non-elaboratable reactants. This is done in the second step of a synthesis route.

    INPUT:
        elab_reactant: mol object of the reactant to elaborate
        product_df: df of the products as produced from step 1.
        ori_reaction: str of the reaction name
        which_elab: int, 1 or 2 corresponding to the reactant label that you want to elaborate.

    OUTPUT:
        resultsdf: dataframe of final products of the elaborations.
    """
    if reaction_name not in REACTIONS_NAMES:
        print(f"Do not have SMARTS for this reaction: {reaction_name}\n "
              f"Please provide the SMARTS.\n")
        return None

    # 1. Make analogs_react df for elaborateable reactant
    no_elab = None
    if which_elab == 1:
        print('ELABORATING REACTANT 1')
        reactant_elab_smarts = REACTANT1_DICT[reaction_name]
        no_elab = 2
        reactant_no_elab_smarts = REACTANT2_DICT[reaction_name]
    else:
        print('ELABORATING REACTANT 2')
        reactant_elab_smarts = REACTANT2_DICT[reaction_name]
        no_elab = 1
        reactant_no_elab_smarts = REACTANT1_DICT[reaction_name]
    analogs_react = processSubMolReactant(elab_reactant, reaction_name, reactant_elab_smarts, reaction_atoms=None,
                                          similarityThr=similarityThr, structuralThr=structuralThr,
                                          reactant_Num=which_elab, resultsDirs=resultsDirs)

    # 2. Get df of already elaborated reactant (that was product of previous step) into same format as analogs_react df
    product_react = productsToReact(product_df, reactant_num=no_elab, output_name=output_name, reactant_smarts=reactant_no_elab_smarts, reaction_name=reaction_name)

    # 3. React back together
    resultsdf = reactBackTogether(original_mol, analogs_react, product_react, reaction_name,
                                  output_name=output_name, max_combinations=max_combinations)

    return resultsdf

def searchReactantAnalogues(originalMol, reactant1, reactant2, resultsDirs, step_num,
                            r1_label=None,
                            r2_label=None,
                            ori_expansionAtomIdx=None,
                            similarityThr=config.SIMILARITY_SEARCH_THR,
                            structuralThr=config.STRUCTURAL_SCORE_THR, struct_score=False, ori_reaction=None,
                            output_name='test', fragmen_info: dict =None):
    """
    Find analogues keeping the atoms involved in the reaction fixed.
    Author: Kate Fieseler

    :param originalMol: original molecule
    :param reactant1:
    :param reactant2:
    :param ori_expansionAtomIdx:
    :param similarityThr:
    :param structuralThr:
    :param struct_score: True if want to use structural score
    :param ori_reaction: original reaction that links both reactants
    :param resultsDirs: list of results directories
    :return:
    """

    print('This is the original molecule:', Chem.MolToSmiles(originalMol))
    print('This is the first reactant:', Chem.MolToSmiles(reactant1))
    print('This is the second reactant:', Chem.MolToSmiles(reactant2))

    # IMPORTANT: ADDING HYDROGENS TO ALL MOLECULES
    Chem.AddHs(originalMol)
    Chem.AddHs(reactant1)
    Chem.AddHs(reactant2)

    # Find attachment points on reactants
    # Returns a list of tuples, (attachmentIdx_whole, attachmentIdx_subMol)
    # Find attachment points on reactants given reaction name
    reactant1_attachmentIdxs = find_attachmentIdxs_fromMolAndSubstruct(originalMol, reactant1, reaction=ori_reaction)
    reactant2_attachmentIdxs = find_attachmentIdxs_fromMolAndSubstruct(originalMol, reactant2, reaction=ori_reaction)

    print('reactant1_attachmentIdxs:', reactant1_attachmentIdxs)
    print('reactant2_attachmentIdxs:', reactant2_attachmentIdxs)

    if len(reactant1_attachmentIdxs) == 0 or len(reactant2_attachmentIdxs) == 0:
        print(f'No attachment points found for reactants... Stopping for this reaction {ori_reaction}')
        return None

    # Find atoms involved in the reaction
    # SMARTS pattern
    if ori_reaction is None:
        print('No reaction provided, need reaction name to find atoms involved in the reaction')
        return NotImplementedError  # TODO: Implement this!
    else:
        print('This is the original reaction:', ori_reaction)
        # Find atoms involved in the reaction
        # TODO: SEND THIS METHOD TO GRAVEYARD
        # react1_atoms = findReactionAtoms(originalMol, reactant1, ori_reaction, resultsDir, 1, reactant1_attachmentIdxs)
        # react2_atoms = findReactionAtoms(originalMol, reactant2, ori_reaction, resultsDir, 2, reactant2_attachmentIdxs)

        (r1_smarts, reactant1, react1_atoms), (r2_smarts, reactant2, react2_atoms) = findReactionAtoms_bothReactants(originalMol, reactant1, reactant2, ori_reaction, resultsDirs, reactant1_attachmentIdxs, reactant2_attachmentIdxs)

    # Perform analogue search of reactants, keeping attachment point fixed
    analogs_react1 = processSubMolReactant(reactant1, ori_reaction, r1_smarts, react1_atoms, similarityThr,
                                           structuralThr, reactant_Num=1,
                                           resultsDirs=resultsDirs)  # are the reacting atoms important here if indexing will change depending on the molecule?
    analogs_react2 = processSubMolReactant(reactant2, ori_reaction, r2_smarts, react2_atoms, similarityThr,
                                           structuralThr, reactant_Num=2, resultsDirs=resultsDirs)

    # Perform reaction search of reactant analogues (including PAINS filter)
    resultsdf = reactBackTogether(originalMol, analogs_react1, analogs_react2, ori_reaction, output_name=output_name, max_combinations=10_000)

    # Add fragmen_info to resultsdf
    if fragmen_info is not None:
        for key, value in fragmen_info.items():
            resultsdf[key] = value

    return resultsdf

