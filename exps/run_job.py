#!/usr/bin/env python3
"""
Run the job with different variable changes.
"""
import sys
import os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
import argparse
from syndirella.cobblers_workshop._base import CobblersWorkshop
from syndirella.slipper._base import Slipper
from typing import (List, Dict, Tuple, Union, Any)

def config_parser():
    parser = argparse.ArgumentParser(description="Run the job with different variable changes.")
    parser.add_argument('--dundee', action='store_true', help="Run pipeline with dundee filtering.")
    parser.add_argument('--kclust', action='store_true', help="Run pipeline with k-means clustering.")
    parser.add_argument('--output', type=str, help="Output directory for the pipeline.")
    parser.add_argument('--template', type=str, help="Absolute path to template for placements.")
    parser.add_argument('--hits', type=str, help="Absolute path to hits for placements. Can be either"
                                                 " .sdf or .mol.")
    return parser

def main():
    parser = config_parser()
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    product = 'c1ccc(-c2cccc3cc[nH]c23)nc1'
    reactants = [('OB(O)c1cccc2cc[nH]c12', 'Brc1ccccn1')]
    reaction_names = ['Sp2-sp2_Suzuki_coupling']
    num_steps = 1
    output_dir: str = '/Users/kate_fieseler/PycharmProjects/syndirella/exps/dundee_noclust' #settings['output']
    filter: bool = True #settings['dundee']
    cluster: bool = False #settings['kclust']
    template = '/Users/kate_fieseler/PycharmProjects/syndirella/exps/fragments/x0310_template.pdb' #settings['template']
    hits = '/Users/kate_fieseler/PycharmProjects/syndirella/exps/fragments/clean_hits.sdf' #settings['hits']
    hits_names = ['x0566_0A']
    batch_num = 3

    # Run the pipeline
    cobblers_workshop = CobblersWorkshop(product, reactants, reaction_names, num_steps,
                                         output_dir, filter)
    final_library = cobblers_workshop.get_final_library()
    final_library.save()
    slipper = Slipper(final_library, template, hits, hits_names, batch_num)
    products = slipper.get_products()
    placements = slipper.place_products()
    slipper.clean_up_placements()

    print('Done!')

if __name__ == '__main__':
    main()