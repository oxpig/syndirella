#!/usr/bin/env python3
"""
syndirella.cli_defaults.py

This file is used to store the default settings for the command line interface.
"""
import os

# Get the scaffold path for the syndirella package from an environment variable
# If the environment variable is not set, a default path is used
# syndirella_base_path = os.environ['SYNDIRELLA_BASE_PATH']
syndirella_base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# Append the specific paths to the scaffold path
print(syndirella_base_path)
rxn_smarts_path = os.path.join(syndirella_base_path, "syndirella/constants/RXN_SMIRKS_CONSTANTS.json")
reactant_filters_path = os.path.join(syndirella_base_path, "syndirella/constants/REACTANT_FILTER_CONSTANTS.json")
additional_rxn_options_path = os.path.join(syndirella_base_path, "syndirella/constants/ADDITIONAL_RXN_OPTIONS.json")

cli_default_settings = dict(
    ranking=os.environ.get('SYNDIRELLA_RANKING', 'num_atom_diff'),
    cutoff=float(os.environ.get('SYNDIRELLA_CUTOFF', 5.0)),
    rxn_smarts_path=rxn_smarts_path,
    reactant_filters_path=reactant_filters_path,
    additional_rxn_options_path=additional_rxn_options_path,
    supressed_exceptions=(Exception,)
)
