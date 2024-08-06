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


def config_parser():
    """
    Configure command-line argument parsing.
    """
    parser = argparse.ArgumentParser(description="Run the Syndirella pipeline with specified configurations.")
    parser.add_argument('--input', type=str, required=True, help="Input .csv file path for the pipeline.")
    parser.add_argument('--output', type=str, required=True, help="Output directory for the pipeline results.")
    parser.add_argument('--templates', type=str, required=True, help="Absolute path to a directory containing the "
                                                                                  "template(s).")
    parser.add_argument('--hits', type=str, required=True, help="Absolute path to hits for placements (.sdf or .mol).")
    parser.add_argument('--products', type=str, help="Absolute path to products for placements.")
    parser.add_argument('--batch_num', type=int, default=10000, help="Batch number for processing.")
    parser.add_argument('--compound_set', action='store_true', help="Include values in compound_set column in the "
                                                                                      "output.")
    parser.add_argument('--manual', action='store_true', help="Use manual routes for processing.")
    return parser


def run_pipeline(settings: Dict[str, Any], pipeline):
    """
    Run the pipeline with the given settings.
    """
    logging.info('Running the pipeline...')

    pipeline.run_pipeline(
        csv_path=settings['input'],
        output_dir=settings['output'],
        template_dir=settings['templates'],
        hits_path=settings['hits'],
        batch_num=settings['batch_num'],
        additional_columns=settings['additional_columns'],
        manual_routes=settings['manual']
    )

    logging.info('Pipeline execution completed successfully.')


def main():
    """
    Main entry point for the CLI.
    """
    setup_logging()

    parser = config_parser()
    args = parser.parse_args()

    syndirella_base_path = os.getenv('SYNDIRELLA_BASE_PATH')
    if not syndirella_base_path:
        logging.error("Environment variable 'SYNDIRELLA_BASE_PATH' is not set.")
        sys.exit(1)

    # Load the pipeline module
    pipeline = load_pipeline_module(syndirella_base_path, 'syndirella/pipeline.py')

    # Convert argparse Namespace to dictionary
    settings = vars(args)
    if settings['compound_set']:
        settings['additional_columns'] = ['compound_set']
    else:
        settings['additional_columns'] = []

    # Run the pipeline
    try:
        run_pipeline(settings, pipeline)
    except Exception as e:
        tb = traceback.format_exc()
        logging.error(f"An error occurred during pipeline execution: {tb}")
        sys.exit(1)


if __name__ == '__main__':
    main()
