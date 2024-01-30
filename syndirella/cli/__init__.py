#!/usr/bin/env python3

from ..version import __version__

# TODO: Finish the docstring
__doc__ = """
Command line interface to Syndirella.
The the first argument to the command is one of the three option:

NB. Please no spaces in the filenames. If spaces are present quote/double-quote the fle path.

Will mainly output data in dataframes, so you can pipe it to csv, json, etc.
""".strip()

from .parser import SyndirellaParser
from .base import SyndirellaParserBase


def main():
    SyndirellaParser.__doc__ = __doc__
    # init configures the parser, calling it executes it
    SyndirellaParser()()

if __name__ == '__main__':
    main()