from .base import SyndirellaParserBase
from .slipper_fitter import SlipperFitterParser
from .slipper_synthesizer import SlipperSynthesizerParser

class SyndirellaParser(SyndirellaParserBase,
                       SlipperFitterParser,
                       SlipperSynthesizerParser):
    pass