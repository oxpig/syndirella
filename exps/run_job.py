#!/usr/bin/env python3
"""
Run the job with different variable changes.
"""
import sys
import os

import pandas as pd

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
import argparse
from syndirella.cobblers_workshop._base import CobblersWorkshop
from syndirella.cobblers_workshop._library import Library
from syndirella.slipper._base import Slipper
from typing import (List, Dict, Tuple, Union, Any)

def config_parser():
    parser = argparse.ArgumentParser(description="Run the job with different variable changes.")
    parser.add_argument('--output', type=str, help="Output directory for the pipeline.")
    parser.add_argument('--template', type=str, help="Absolute path to template for placements.")
    parser.add_argument('--hits', type=str, help="Absolute path to hits for placements. Can be either"
                                                 " .sdf or .mol.")
    parser.add_argument('--products', type=str, help="Absolute path to products for placements.")
    return parser

def main():
    parser = config_parser()
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    product = 'c1ccc(-c2cccc3cc[nH]c23)nc1'
    reactants = [('OB(O)c1cccc2cc[nH]c12', 'Brc1ccccn1')]
    reaction_names = ['Sp2-sp2_Suzuki_coupling']
    num_steps = 1
    output_dir: str = settings['output']
    template = settings['template']
    hits = settings['hits']
    hits_names = ['x0566_0A']
    batch_num = 10
    final_products_library_csv_path = settings['products']
    final_products = pd.read_csv(final_products_library_csv_path, index_col=0)

    # Run placement
    final_library = Library.load(output_dir)
    slipper = Slipper(final_library, template, hits, hits_names, batch_num, False)
    slipper.final_products_csv_path = final_products_library_csv_path
    slipper.products = final_products
    placements = slipper.place_products()
    slipper.clean_up_placements()

    print('Done!')

if __name__ == '__main__':
    main()