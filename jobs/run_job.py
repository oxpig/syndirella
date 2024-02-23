#!/usr/bin/env python3
"""
Run the pipeline 🏁
"""
import sys
import os

import pandas as pd

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
import argparse
from typing import (List, Dict, Tuple, Union, Any)

from syndirella.pipeline import run_pipeline

def config_parser():
    parser = argparse.ArgumentParser(description="Run the job with different variable changes.")
    parser.add_argument('--input', type=str, help="Input .csv for the pipeline.")
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
    input_csv_path: str = settings['input']
    output_dir: str = settings['output']
    template = settings['template']
    hits = settings['hits']
    batch_num = 100
    additional_info = ['compound_set']
    manual_routes = True

    # Run pipeline
    print('Running pipeline...')
    run_pipeline(input_csv_path,
                 output_dir,
                 template,
                 hits,
                 batch_num,
                 additional_info,
                 manual_routes)
    print('Done!')

if __name__ == '__main__':
    main()