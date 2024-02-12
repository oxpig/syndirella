#!/usr/bin/env python3
"""
syndirella._cli_defaults.py

This file is used to store the default settings for the command line interface.
"""
import os

cli_default_settings = dict(
    ranking=os.environ.get('SYNDIRELLA_RANKING', 'num_atom_diff'),
    cutoff=float(os.environ.get('SYNDIRELLA_CUTOFF', 5.0)),
    rxn_smarts_path=os.environ.get('SYNDIRELLA_RXN_SMARTS_PATH',
                                   "/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/"
                                   "constants/RXN_SMARTS_CONSTANTS.json"),
    reactant_filters_path=os.environ.get('SYNDIRELLA_REACTANT_FILTER_PATH',
                                         "/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/"
                                         "constants/REACTANT_FILTER_CONSTANTS_V2.json"),
    supressed_exceptions=(Exception,),
)
