#!/usr/bin/env python3
"""
CLI script to run the Syndirella pipeline or justretroquery 🏁
"""
import argparse
import cProfile
import importlib.util
import logging
import os
import sys
import traceback
from typing import Dict, Any

from syndirella.constants import DatabaseSearchTool, RetrosynthesisTool, DEFAULT_DATABASE_SEARCH_TOOL, DEFAULT_RETROSYNTHESIS_TOOL


class SyndirellaHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Custom help formatter that shows command arguments more clearly."""
    
    def _format_action(self, action):
        # Get the help from the parent class
        help_text = super()._format_action(action)
        
        # If this is a subparser action, format it specially
        if isinstance(action, argparse._SubParsersAction):
            # Get the subparsers
            subparsers = []
            for choice, subparser in action.choices.items():
                # Get the description if available
                description = getattr(subparser, 'description', '')
                subparsers.append(f"  {choice:<15} {subparser.description or subparser.help}")
            
            # Format the subparsers
            help_text = f"{action.help}:\n" + "\n".join(subparsers)
        
        return help_text


def setup_logging(level=logging.INFO):
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )


def load_pipeline_module(syndirella_base_path: str, module_relative_path: str):
    pipeline_module_path = os.path.join(syndirella_base_path, module_relative_path)
    spec = importlib.util.spec_from_file_location("syndirella.pipeline", pipeline_module_path)
    pipeline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)
    return pipeline


def validate_api_credentials(settings: Dict[str, Any]):
    retro_tool = RetrosynthesisTool.from_string(settings.get('retro_tool', DEFAULT_RETROSYNTHESIS_TOOL.value))
    db_search_tool = DatabaseSearchTool.from_string(settings.get('db_search_tool', DEFAULT_DATABASE_SEARCH_TOOL.value))
    
    if retro_tool == RetrosynthesisTool.MANIFOLD:
        if not os.environ.get('MANIFOLD_API_KEY'):
            logging.error("MANIFOLD_API_KEY environment variable not set.")
            sys.exit(1)
        if not os.environ.get('MANIFOLD_API_URL'):
            logging.error("MANIFOLD_API_URL environment variable not set.")
            sys.exit(1)
    
    if db_search_tool == DatabaseSearchTool.ARTHOR:
        if not os.environ.get('ARTHOR_API_URL'):
            logging.warning("ARTHOR_API_URL environment variable not set. Using default: https://arthor.docking.org")


def config_parser(syndirella_base_path: str):
    parser = argparse.ArgumentParser(prog="syndirella",
                                     description="Run the Syndirella pipeline with specified configurations.",
                                     epilog=f"Syndirella is installed at {syndirella_base_path} \n",
                                     formatter_class=SyndirellaHelpFormatter)
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add setup-aizynth command
    setup_parser = subparsers.add_parser('setup-aizynth', 
                                        help='Setup AiZynthFinder data and configuration',
                                        description='Automatically download AiZynthFinder data and create configuration file.')
    
    pipeline_parser = subparsers.add_parser('run', 
                                           help='Run the main Syndirella pipeline',
                                           description='Run the full Syndirella pipeline with specified input files and parameters.')
    pipeline_parser.add_argument('-i', '--input', type=str, required=True, help="Input .csv file path for the pipeline.")
    pipeline_parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for the pipeline results.")
    pipeline_parser.add_argument('-t', '--templates', type=str, required=False,
                        help="Absolute path to a directory containing the template(s).")
    pipeline_parser.add_argument('--hits_path', type=str, required=False,
                        help="Absolute path to hits_path for placements (.sdf or .mol).")
    pipeline_parser.add_argument('--products', type=str, required=False, help="Absolute path to products for placements.")
    pipeline_parser.add_argument('--batch_num', type=int, default=10000, help="Batch number for processing.")
    pipeline_parser.add_argument('--manual', action='store_true', help="Use manual routes for processing.")
    pipeline_parser.add_argument('--only_scaffold_place', action='store_true',
                        help="Only place scaffolds. Do not continue to elaborate.")
    pipeline_parser.add_argument('--scaffold_place_num', type=int, default=5,
                        help="Number of times to attempt scaffold placement.")
    pipeline_parser.add_argument('--retro_tool', type=str, default=DEFAULT_RETROSYNTHESIS_TOOL.value,
                        choices=[tool.value for tool in RetrosynthesisTool],
                        help="Retrosynthesis tool to use.")
    pipeline_parser.add_argument('--db_search_tool', type=str, default=DEFAULT_DATABASE_SEARCH_TOOL.value,
                        choices=[tool.value for tool in DatabaseSearchTool],
                        help="Database search tool to use.")
    pipeline_parser.add_argument('--profile', action='store_true', help="Run the pipeline with profiling.")
    pipeline_parser.add_argument('--atom_diff_min', type=int, default=0,
                        help="Minimum atom difference between elaborations and scaffold to keep.")
    pipeline_parser.add_argument('--atom_diff_max', type=int, default=10,
                        help="Maximum atom difference between elaborations and scaffold to keep.")
    pipeline_parser.add_argument('--just_retro', action='store_true', help="Only run retrosynthesis querying of scaffolds.")
    pipeline_parser.add_argument('--no_scaffold_place', action='store_true',
                        help="Do not place scaffolds initially before elaborating.")
    pipeline_parser.add_argument('--elab_single_reactant', action='store_true',
                        help="Only elaborate one reactant per elaboration series.")
    pipeline_parser.add_argument('--reference_db', type=str,
                        help="Path to reference HIPPO database file for superstructure search, must set --db_search_tool to 'hippo'.")
    pipeline_parser.add_argument('--no_assert_scaffold_intra_geom_flatness', action='store_true',
                        help="Don't check scaffold for intra geometry or flatness.")
    pipeline_parser.add_argument('--substructure_check_products', action='store_true',
                                 help="Only keep products that are superstructures of their scaffolds.")
    
    add_reaction_parser = subparsers.add_parser('add-reaction', 
                                               help='Add a new reaction to the library',
                                               description='Add a new reaction SMIRKS to the reaction library with optional parent finding.')
    add_reaction_parser.add_argument('--name', type=str, required=True, help="Name of the new reaction.")
    add_reaction_parser.add_argument('--smirks', type=str, required=True, help="SMIRKS string for the reaction.")
    add_reaction_parser.add_argument('--find_parent', action='store_true', 
                        help="If True, treat as a child reaction and find parent based on similarity.")
    add_reaction_parser.add_argument('--fp_type', type=str, default='maccs_rxn_fp',
                        choices=['maccs_rxn_fp', 'morgan_rxn_fp'],
                        help="Fingerprint type for similarity calculation.")
    add_reaction_parser.add_argument('--threshold', type=float, default=0.2,
                        help="Similarity threshold for finding parent reaction.")
    add_reaction_parser.add_argument('--similarity_metric', type=str, default='jaccard',
                        choices=['jaccard', 'cosine'],
                        help="Similarity metric for finding parent reaction.")
    
    return parser


def show_command_help(command_name: str):
    """Show detailed help for a specific command."""
    syndirella_base_path = os.path.dirname(importlib.util.find_spec('syndirella').origin)
    parser = config_parser(syndirella_base_path)
    
    # Get the subparser for the command
    subparsers = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    if subparsers:
        subparser_action = subparsers[0]
        if command_name in subparser_action.choices:
            subparser_action.choices[command_name].print_help()
        else:
            print(f"Unknown command: {command_name}")
            print("Available commands:")
            for choice in subparser_action.choices:
                print(f"  {choice}")
    else:
        print("No subcommands found")


def setup_aizynthfinder():
    """Setup AiZynthFinder data and configuration."""
    try:
        from syndirella.aizynth.AiZynthManager import AiZynthManager
        
        # Create a temporary manager to trigger setup
        manager = AiZynthManager(auto_setup=True)
        config_path = manager.config_file
        logging.info(f"AiZynthFinder setup complete. Config: {config_path}")
        logging.info("Note: This environment variable is set for the current session only.")
        logging.info("To make it permanent, add this line to your shell profile (~/.bashrc, ~/.zshrc, etc.):")
        logging.info(f'export AIZYNTH_CONFIG_FILE="{config_path}"')

        
    except Exception as e:
        logging.error(f"Error setting up AiZynthFinder: {e}")
        logging.error("Please install AiZynthFinder and run 'download_public_data' manually.")
        sys.exit(1)


def run_pipeline(settings: Dict[str, Any], pipeline):
    logging.info('Running the pipeline...')
    pipeline.run_pipeline(settings)
    logging.info('Pipeline execution completed successfully.')


def run_justretroquery(settings: Dict[str, Any], justretroquery):
    logging.info('Running justretroquery...')
    justretroquery.run_justretroquery(settings)
    logging.info('Justretroquery execution completed successfully.')


def main():
    setup_logging()
    
    # Find syndirella package location more robustly
    try:
        # First try the standard approach
        spec = importlib.util.find_spec('syndirella')
        if spec and spec.origin:
            syndirella_base_path = os.path.dirname(spec.origin)
        else:
            # Fallback: use the package's __file__ attribute
            import syndirella
            syndirella_base_path = os.path.dirname(syndirella.__file__)
    except Exception:
        # Final fallback: assume we're in the project directory
        # This works when running from the project root
        current_dir = os.path.dirname(os.path.abspath(__file__))
        if os.path.basename(current_dir) == 'syndirella':
            syndirella_base_path = current_dir
        else:
            # Look for syndirella directory in the current working directory
            cwd = os.getcwd()
            potential_path = os.path.join(cwd, 'syndirella')
            if os.path.exists(potential_path):
                syndirella_base_path = potential_path
            else:
                raise RuntimeError("Could not determine syndirella package location")

    # Check if user wants help for a specific command
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help'] and len(sys.argv) > 2:
        show_command_help(sys.argv[2])
        return

    parser = config_parser(syndirella_base_path)
    args = parser.parse_args()

    # Handle setup-aizynth command
    if hasattr(args, 'command') and args.command == 'setup-aizynth':
        setup_aizynthfinder()
        return

    if args.command == 'add-reaction':
        try:
            from syndirella.route.SmirksLibraryManager import SmirksLibraryManager
            
            constants_dir = os.path.join(syndirella_base_path, 'constants')
            smirks_path = os.path.join(constants_dir, 'RXN_SMIRKS_CONSTANTS.json')
            uspto_path = os.path.join(constants_dir, 'uspto_template_lookup.json')
            
            manager = SmirksLibraryManager(smirks_path, uspto_path)
            
            result = manager.add_reaction_to_library(
                name=args.name,
                smirks=args.smirks,
                find_parent=args.find_parent,
                fp_type=args.fp_type,
                threshold=args.threshold,
                similarity_metric=args.similarity_metric
            )
            logging.info(f"Successfully added reaction: {result}")
            sys.exit(0)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f"An error occurred while adding reaction: {tb}")
            sys.exit(1)

    if not hasattr(args, 'command') or args.command is None:
        settings = vars(args)
        validate_api_credentials(settings)

        if settings['just_retro']:
            justretroquery = load_pipeline_module(syndirella_base_path, 'justretroquery.py')
            try:
                run_justretroquery(settings, justretroquery)
                sys.exit(0)
            except Exception as e:
                tb = traceback.format_exc()
                logging.error(f"An error occurred during justretroquery execution: {tb}")
                sys.exit(1)

        pipeline = load_pipeline_module(syndirella_base_path, 'pipeline.py')

        if settings['profile']:
            try:
                profiler = cProfile.Profile()
                profiler.enable()
                run_pipeline(settings, pipeline)
                profiler.disable()
                profiler.print_stats(sort='time')
                print('\n\n')
                profiler.print_stats(sort='cumtime')
                sys.exit(0)
            except Exception as e:
                tb = traceback.format_exc()
                logging.error(f"An error occurred during pipeline execution: {tb}")
                sys.exit(1)

        try:
            run_pipeline(settings, pipeline)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f"An error occurred during pipeline execution: {tb}")
            sys.exit(1)

    elif args.command == 'run':
        settings = vars(args)
        validate_api_credentials(settings)

        if settings['just_retro']:
            justretroquery = load_pipeline_module(syndirella_base_path, 'justretroquery.py')
            try:
                run_justretroquery(settings, justretroquery)
                sys.exit(0)
            except Exception as e:
                tb = traceback.format_exc()
                logging.error(f"An error occurred during justretroquery execution: {tb}")
                sys.exit(1)

        pipeline = load_pipeline_module(syndirella_base_path, 'pipeline.py')

        if settings['profile']:
            try:
                profiler = cProfile.Profile()
                profiler.enable()
                run_pipeline(settings, pipeline)
                profiler.disable()
                profiler.print_stats(sort='time')
                print('\n\n')
                profiler.print_stats(sort='cumtime')
                sys.exit(0)
            except Exception as e:
                tb = traceback.format_exc()
                logging.error(f"An error occurred during pipeline execution: {tb}")
                sys.exit(1)

        try:
            run_pipeline(settings, pipeline)
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(f"An error occurred during pipeline execution: {tb}")
            sys.exit(1)


if __name__ == '__main__':
    main()
