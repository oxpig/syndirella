import argparse, os, json, itertools, string
import multiprocessing
from functools import partial

from rdkit import Chem
from .._cli_defaults import cli_default_settings
from .base import set_verbose
from ..SMARTSHandler import SMARTSHandler
from typing import List, Any, Dict
import pandas as pd
from rdkit.Chem import PandasTools
import ast

class SyndirellaParserPipeline:
    def _define_pipeline(self, parser: argparse.ArgumentParser):
        """
        syndirella pipeline -h

        :param parser:
        :return:
        """
        parser.add_argument('-i', '--input', type=str,
                        help=('Path to the input CSV file. The expected CSV structure is:\n'
                              'product (str) - SMILES of product\n'
                              'output_dir (list) - Name of the directories to save results. Usually target id.\n'
                              'num_steps (int) - Number of steps in the route\n'
                              'rxn_order_first_to_last (list(str)) - Reaction name to produce product\n'
                              'reactants (list(tuple)) - Reactants listed in tuples\n'
                              '...\n'))
        parser.add_argument('-o', '--output', help='Output folder',
                            default=cli_default_settings['output'])
        parser.add_argument('-b', '--batch_size', help='Batch size')
        parser.add_argument('-n', '--n_cores', help='Number of cores',
                            default=cli_default_settings['n_cores'],
                            type=int)

        parser.set_defaults(func=self.pipeline)
        return parser

    def convert_to_literal(self, cell: pd.Cell, expected_type):
        try:
            # Safely evaluate the string to the expected literal type
            return ast.literal_eval(cell) if pd.notnull(cell) else cell
        except (ValueError, SyntaxError):
            # Return a default value if conversion fails
            return expected_type()

    def validate_and_convert_columns(self, batch_data: pd.DataFrame, column_specs: Dict):
        for column, expected_type in column_specs.items():
            if expected_type in [list, tuple]:
                # Convert string representation of lists or tuples to actual lists or tuples
                batch_data[column] = batch_data[column].apply(lambda cell: self.convert_to_literal(cell, expected_type))
                # Check if each element in the column is now of the expected type
                if not all(isinstance(cell, expected_type) for cell in batch_data[column] if pd.notnull(cell)):
                    raise TypeError(
                        f"Column '{column}' contains elements that are not of type {expected_type.__name__}.")
            else:
                # For other types like 'int' or 'str', ensure the column is of the expected type
                if not pd.api.types.is_dtype_equal(batch_data[column].dtype, 'object'):
                    batch_data[column] = batch_data[column].astype(expected_type)
                else:
                    # If column is object type, check all values are of the expected Python type
                    if not all(isinstance(cell, expected_type) for cell in batch_data[column] if pd.notnull(cell)):
                        raise TypeError(f"Column '{column}' contains elements that are not of type {expected_type}.")

    def process_example(self, example: pd.Series, shared_resources):
        """
        ENTRY POINT TO REST OF SYNDIRELLA.
        """
        smarts_handler = SMARTSHandler()
        cobbler_bench = CobblerBench(product=example['product'], output_dir=example['output_dir'], num_steps=example['num_steps'],
                                     reaction_names=example['rxn_order_first_to_last'], reactants=example['reactants'],
                                     smarts_handler=smarts_handler)
        '''
            GOAL OF PROCESS_EXAMPLE:
            def process_example(self, example: pd.Series, shared_resources):
                cobbler_workshop = CobblerWorkshop() # will automatically know if single or multi
                library = cobbler_workshop.get_final_library() # will return the final product library
                reaction = library.reaction() # will be the final reaction object

                first_library = cobbler_library.get_first_library() # will return the first library
                first_reaction = first_library.reaction() # will be the first reaction object

                slipper = Slipper(library) # will be the final slipper object, requires cobbler_library to be instantiated
                slipper.get_products() # will return the final slippers using SlipperSynthesizer class
                slipper.place_products() # will return the final placed slippers using SlipperFitter class
        '''
        pass

    def process_batch(self, batch_data: pd.DataFrame, settings: Dict[str, Any], batch_index: int, shared_resources):
        # Define required columns and their expected types
        column_specs = {
            'product': str,
            'output_dir': list,
            'num_steps': int,
            'rxn_order_first_to_last': list,
            'reactants': list  # Assuming that you have a way to verify that the list contains tuples
        }
        # Validate if required columns are present and convert types where necessary
        self.validate_and_convert_columns(batch_data, column_specs)
        # If validation and conversion are successful, process each example in the batch
        for index, example in batch_data.iterrows():
            self.process_example(example, shared_resources)
        print(f"Batch {batch_index} processed.")

    def pipeline(self, args: argparse.Namespace):
        """
        Performs the whole syndirella pipeline.
        Identifies which reactant is which in the reaction SMARTS of the reaction template.
        Gets analogues of those reactnts using Library class.
        Combines the analogues to get the products using the SlipperSynthesizer class
        Places the products in the template using the SlipperFitter class
        """
        set_verbose(args.verbose)
        settings: Dict[str, Any] = vars(args)
        # Read in input csv
        df: pd.DataFrame = pd.read_csv(settings['input'], sep=',', header=0)
        # Define batch size
        batch_size = settings['batch_size']  # Ensure 'batch_size' is defined in args
        # Segment dataframe into batches
        batches = [df[i:i + batch_size] for i in range(0, df.shape[0], batch_size)]
        # Process each batch
        if args.n_cores > 1:
            # Use multiprocessing to process batches in parallel if more than one core is specified
            pool = multiprocessing.Pool(args.n_cores)
            # Partial function to pass the additional args
            func = partial(self.process_batch, settings=settings)
            pool.starmap(func, enumerate(batches))
            pool.close()
            pool.join()
        else:
            # Process batches sequentially
            for index, batch_data in enumerate(batches):
                self.process_batch(batch_data, args, index)
        print("All batches processed.")

