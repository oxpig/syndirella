#!/usr/bin/env python3
"""
CLI script to run the Syndirella pipeline 🏁
"""
import sys
import os
import argparse
import importlib.util
import logging
from typing import Dict, Any
import traceback
import cProfile

def setup_logging(level=logging.INFO):
    """
    Set up logging configuration
    """
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )


def load_pipeline_module(syndirella_base_path: str, module_relative_path: str):
    """
    Dynamically load the pipeline module from a specified path.
    """
    pipeline_module_path = os.path.join(syndirella_base_path, module_relative_path)
    spec = importlib.util.spec_from_file_location("syndirella.pipeline", pipeline_module_path)
    pipeline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)
    return pipeline


def config_parser(syndirella_base_path: str):
    """
    Configure command-line argument parsing.
    """
    parser = argparse.ArgumentParser(prog="syndirella",
                                     description="Run the Syndirella pipeline with specified configurations.",
                                     epilog=f"Syndirella is installed at {syndirella_base_path} \n")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input .csv file path for the pipeline.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for the pipeline results.")
    parser.add_argument('-t', '--templates', type=str, required=True, help="Absolute path to a directory containing the "
                                                                                   "template(s).")
    parser.add_argument('--hits_path', type=str, required=True, help="Absolute path to hits_path for placements (.sdf or .mol).")
    parser.add_argument('--metadata', type=str, required=True, help="Absolute path to metadata for placements.")
    parser.add_argument('--products', type=str, required=False, help="Absolute path to products for placements.")
    parser.add_argument('--batch_num', type=int, default=10000, help="Batch number for processing.")
    parser.add_argument('--manual', action='store_true', help="Use manual routes for processing.")
    parser.add_argument('--profile', action='store_true', help="Run the pipeline with profiling.")

    return parser


def run_pipeline(settings: Dict[str, Any], pipeline):
    """
    Run the pipeline with the given settings.
    """
    logging.info('Running the pipeline...')

    pipeline.run_pipeline(csv_path=settings['input'],
                          output_dir=settings['output'],
                          template_dir=settings['templates'],
                          hits_path=settings['hits_path'],
                          metadata_path=settings['metadata'],
                          batch_num=settings['batch_num'],
                          additional_columns=['compound_set'], # Will always be compound_set
                          manual_routes=settings['manual'])

    logging.info('Pipeline execution completed successfully.')


def main():
    """
    Main entry point for the CLI.
    """
    setup_logging()
    syndirella_base_path = os.path.dirname(importlib.util.find_spec('syndirella').origin)

    parser = config_parser(syndirella_base_path)
    args = parser.parse_args()

    # check for MANIFOLD API key
    if not os.environ.get('MANIFOLD_API_KEY'):
        logging.error("MANIFOLD_API_KEY environment variable not set.")
        sys.exit(1)

    # Load the pipeline module
    pipeline = load_pipeline_module(syndirella_base_path, 'pipeline.py')

    # Convert argparse Namespace to dictionary
    settings = vars(args)

    if settings['profile']:
        profiler = cProfile.Profile()
        profiler.enable()
        run_pipeline(settings, pipeline)
        profiler.disable()
        profiler.print_stats(sort='time')
        print('\n\n')
        profiler.print_stats(sort='cumtime')
    else:
        # Run the pipeline
        try:
            run_pipeline(settings, pipeline)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f"An error occurred during pipeline execution: {tb}")
            sys.exit(1)


if __name__ == '__main__':
    main()
