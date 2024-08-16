#!/usr/bin/env python3
"""
error.py

Possible errors when using syndirella.
"""

class NoReactants(Exception):
    pass

class NoSynthesisRoute(Exception):
    pass

class ReactionError(Exception):
    pass

class ScaffoldPlacementError(Exception):
    def __init__(self,
                 scaffold: str):
        self.scaffold = scaffold

    def __str__(self):
        return f'{self.scaffold} could not be successfully placed.'