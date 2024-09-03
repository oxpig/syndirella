#!/usr/bin/env python3
"""
Goal: Profile on 3 different scaffold compounds placing 10 elaborations each to get good statistics on the pipeline.

"""
import argparse
import cProfile
import sys
import os
import importlib.util

def config_parser():
    parser = argparse.ArgumentParser(description="Run the job with different variable changes.")
    parser.add_argument('--input', type=str, help="Input .csv name.")
    parser.add_argument('--output', type=str, help="Output name of file with profiling information (Does"
                                                   "not need to include path, just name of file you'd like).")
    parser.add_argument('--base_dir', type=str, help="Base directory for the pipeline.")
    return parser

def main():
    parser = config_parser()
    # load
    settings = vars(parser.parse_args())
    input_csv_name = settings['input']
    output_file = settings['output']
    base_dir = settings['base_dir']

    # Set the paths
    input_csv_path = os.path.join(base_dir, input_csv_name)
    output_dir = os.path.join(base_dir, 'outputs')
    template = os.path.join(base_dir, 'Ax0310_relaxed_apo.pdb')
    hits = os.path.join(base_dir, 'A71EV2A_combined.sdf')
    batch_num = 100
    additional_info = ['compound_set']
    manual_routes = True

    # Set output file name
    output_file = os.path.join(base_dir, output_file)

    # Syndirella path
    sys.path.append(os.environ['SYNDIRELLA_BASE_PATH'])
    # Get the scaffold path and construct the full path to the module
    syndirella_path = os.environ['SYNDIRELLA_BASE_PATH']
    pipeline_module_path = os.path.join(syndirella_path, 'syndirella/pipeline.py')
    # Load the module specified by the path
    spec = importlib.util.spec_from_file_location("syndirella.pipeline", pipeline_module_path)
    pipeline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)

    def run_pipeline():
        # Check which shell
        print('The shell being used is: ', os.getenv('SHELL'))
        # Run pipeline
        print('Running pipeline...')
        pipeline.run_pipeline(csv_path=input_csv_path,
                              output_dir=output_dir,
                              template_path=template,
                              hits_path=hits,
                              batch_num=batch_num,
                              additional_columns=additional_info,
                              manual_routes=manual_routes)
        print('Done!')

    # Redirect stdout to a file
    # original_stdout = sys.stdout
    # with open(output_file, 'w') as f:
    #     sys.stdout = f
    profiler = cProfile.Profile()
    profiler.enable()
    run_pipeline()
    profiler.disable()
    profiler.print_stats(sort='time')
    print('\n\n')
    profiler.print_stats(sort='cumtime')

    # Reset stdout to original
    #sys.stdout = original_stdout

if __name__ == '__main__':
    main()
