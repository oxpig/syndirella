import argparse, logging
from typing import Optional, List


class SyndirellaParserBase:
    """
    See main module docstring for usage.
    """

    subcommand_help = '''Actions: slipperfitter, slippersynthesizer'''
    slipperfitter_help = '''slipperfitter'''
    slippersynthesizer_help = '''slippersynthesizer'''
    pipeline_help = '''pipeline'''
    utils_help = '''utils'''

    def __init__(self):
        """
        Defines the parser.
        Calling executes it.
        """
        self.parser = argparse.ArgumentParser(description=self.__doc__)
        subparsers = self.parser.add_subparsers(title='subcommands', help=self.subcommand_help)
        slippersynthesizer_parser = subparsers.add_parser('slippersynthesizer', help=self.slippersynthesizer_help)
        self._define_slippersynthesizer(slippersynthesizer_parser)
        slipperfitter_parser = subparsers.add_parser('slipperfitter', help=self.slipperfitter_help)
        self._define_slipperfitter(slipperfitter_parser)
        pipeline_parser = subparsers.add_parser('pipeline', help=self.pipeline_help)
        self._define_pipeline(pipeline_parser)
        utils_parser = subparsers.add_parser('utils', help=self.utils_help)
        self._define_utils(utils_parser)

    def __call__(self, cli_override: Optional[List[str]] = None):
        if cli_override:
            args = self.parser.parse_args(cli_override)
        else:
            args = self.parser.parse_args()
        args.func(args)

    def _define_slippersynthesizer(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')

    def _define_slipperfitter(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')

    def _define_utils(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')

    def _define_pipeline(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')


def set_verbose(count):
    log_map = {None: logging.FATAL,
               0: logging.CRITICAL,
               1: logging.ERROR,
               2: logging.WARNING,
               3: logging.INFO,
               4: logging.DEBUG}
    #Victor.enable_stdout(log_map.get(count, logging.DEBUG))
    #Victor.capture_logs()

"""
def _add_common(parser: argparse.ArgumentParser, **settings):
    if settings.get('verbose', True):
        parser.add_argument('-v', '--verbose', action="count", help='verbose')
    if settings.get('hits_path', False):
        parser.add_argument('-i', '--hits_path', nargs='+', required=True, help='hit mol files')
    if settings.get('output', False):
        parser.add_argument('-o', '--output', default='.', help='output root folder')
    if settings.get('name', False):
        parser.add_argument('-n', '--name', default='fragmenstein', help='output name of molecule')
    if settings.get('template', False):
        parser.add_argument('-t', '--template', required=True, help='template PDB file')
"""