import argparse
import pandas as pd
from .base import SyndirellaParserBase, _add_common, set_verbose

class SlipperSynthesizerParser:
    def _define_slippersynthesizer(self, parser: argparse.ArgumentParser):
        """
        syndirella slippersynthesizer -h

        :param parser:
        :return:
        """
        parser.add_argument('-i', '--input', help='Input file', required=True)
        parser.add_argument('-o', '--output', help='Output file', required=True)
        parser.add_argument('-v', '--verbose', help='Verbosity level', default=0, type=int)
        parser.add_argument('-s', '--smiles', help='Smiles column name', default='SMILES')
        parser.add_argument('-t', '--template', help='Template column name', default='Template')
        parser.add_argument('-r', '--reaction', help='Reaction column name', default='Reaction')
        parser.add_argument('-c', '--cutoff', help='Cutoff', default=0.5, type=float)
        parser.add_argument('-m', '--max', help='Max number of products', default=5, type=int)
        parser.add_argument('-n', '--n_cores', help='Number of cores', default=1, type=int)
        parser.add_argument('-p', '--parallel', help='Parallel', default=False, type=bool)
        parser.add_argument('-d', '--database', help='Database', default='synth')
        parser.add_argument('-a', '--append', help='Append', default=False, type=bool)
        parser.add_argument('-f', '--force', help='Force', default=False, type=bool)
        parser.add_argument('-b', '--batch_size', help='Batch size', default=100, type=int)
        parser.add_argument('-e', '--error', help='Error file', default='error.csv')
        parser.add_argument('-l', '--log', help='Log file', default='log.csv')
        parser.add_argument('-x', '--index', help='Index column name', default='Index')
        parser.add_argument('-z', '--zip', help='Zip', default=False, type=bool)
        parser.set_defaults(func=self._slippersynthesizer)