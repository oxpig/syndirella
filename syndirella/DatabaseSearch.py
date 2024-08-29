#!/usr/bin/env python3
"""
DatabaseSearch.py

This module contains the functionality for a database search. Can do it via Postera or another database (smallworld, etc.).
"""
from syndirella.route.Reaction import Reaction
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple, Optional)
from rdkit import Chem
from abc import ABC, abstractmethod

class DatabaseSearch(ABC):
    """
    This class contains information about the database search. It will perform the database search using the
    database_search_function. It will also store the results of the database search as a .csv file.
    """
    def __init__(self):
        self.url = None
        self.api_key = None

    @abstractmethod
    def perform_database_search(self,
                                reactant: Chem.Mol,
                                reaction_name: str,
                                search_type: str):
        """
        This function is used to perform the database search using the database_search_function.
        """
        pass

    @abstractmethod
    def get_resp_json(self, url: str, data: Dict = None) -> Optional[Dict]:
        """
        Directly get the response json from a request.
        """
        pass

    @abstractmethod
    def save_database_search(self):
        """
        This function is used to save the results of the database search as a .csv file.
        """
        pass

    @abstractmethod
    def load_database_search(self):
        """
        This function is used to load the results of the database search from a .csv file.
        """
        pass
