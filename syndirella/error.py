#!/usr/bin/env python3
"""
error.py

Possible errors when using syndirella.
"""

class ReactionError(Exception):
    pass

class RouteError(Exception):
    def __init__(self,
                 base: str):
        self.base = base

    def __str__(self):
        return f'{self.base} could not be successfully placed.'

