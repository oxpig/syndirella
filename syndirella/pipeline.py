#!/usr/bin/env python3
"""
syndirella.pipeline.py

This module contains the Pipeline class. This is the main class that performs the whole elaboration process.
"""
import os
from typing import List, Dict, Tuple, Union, Optional

class Pipeline:
    """
    The Pipeline class is used to perform the whole elaboration process.
    """
    def __init__(self,
             csv_path: str):
        self.csv_path: str = self.assert_csv_path(csv_path)

    def assert_csv_path(self, csv_path: str) -> str:

