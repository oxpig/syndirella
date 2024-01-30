#!/usr/bin/env python3
"""
_postera.py

This module contains the functionality for a Postera search.
"""
from ._database_search import DatabaseSearch
from ._reaction import Reaction
import os
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict, Tuple, Optional)
from rdkit import Chem
import requests
import json
from collections import OrderedDict
from rdkit import DataStructs
from rdkit.Chem import AllChem

class Postera(DatabaseSearch):
    """
    This class contains information about the Postera search. It will perform the Postera search using the
    perform_database_search function. It will also store the results of the Postera search as a .csv file.
    """
    def __init__(self, reaction: Reaction, output_dir: str, search_type: str):
        super().__init__(reaction, output_dir, search_type)
        self.url = "https://api.postera.ai"
        self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.search_type = search_type

    def perform_database_search(self, reactant: Chem.Mol):
        """
        This function is used to perform the Postera search using the database_search_function.
        """
        smiles = Chem.MolToSmiles(reactant)
        if self.search_type == "superstructure":
            hits: List[str] = self.perform_superstructure_search(smiles)
        elif self.search_type == "ring_replace":
            hits: List[str] = self.perform_ring_replace_search(smiles)
        else:
            raise ValueError("Error, library_type must be either 'superstructure' or 'ring_replace'.")
        filtered_hits: List[str] = self.filter_out_hits(hits, reactant)
        print(f'Found {len(filtered_hits)} hits_path for {smiles}')
        return hits

    def filter_out_hits(self, hits: List[str], reactant: Chem.Mol) -> List[str]:
        """
        This function is used to filter out hits_path that have
            chirality specification
            repeat hits_path
            non-abundant isotopes
            original reactant
        """
        hits_mols = [Chem.MolFromSmiles(hit) for hit in hits]
        filtered_mols = list(OrderedDict((DataStructs.BitVectToText(AllChem.GetMorganFingerprintAsBitVect(
            mol, 2)), mol) for mol in hits_mols).values())
        filtered_hits = [Chem.MolToSmiles(mol) for mol in filtered_mols]
        filtered_hits = [hit for hit in filtered_hits if "@" not in hit] # remove chirality specification
        filtered_hits = [hit for hit in filtered_hits if not ('15' in hit) or ('13' in hit)] # TODO: Test if this works
        # Remove original reactant if is in the list of hits_path
        r_smiles = Chem.MolToSmiles(reactant)
        if r_smiles in filtered_hits: filtered_hits.remove(r_smiles)
        return filtered_hits

    def perform_superstructure_search(self, smiles: str, max_pages: int = 10) -> List[str]:
        """
        This function is used to perform the Postera superstructure search.
        """
        print(f"Running superstructure search for {smiles}")
        superstructure_hits = self.get_search_results(
            url=f'{self.url}/api/v1/superstructure/',
            data={
                'smiles': smiles,
                # You could pass in more filters here
            },
            max_pages=max_pages,
        )
        return [r["smiles"] for r in superstructure_hits]

    def get_search_results(self, url: str, data: Dict[str, Any],
                           max_pages: int = 5, page: int = 1) -> List[Dict]:
        """
        Recursively get all pages for the endpoint until reach null next page or
        the max_pages threshold.
        """
        # List where we will gather up all the hits_path.
        all_hits = []
        data = {
            **data,
            'page': page,
        }
        response = self.get_resp_json(url, data)
        if response is None:
            return all_hits
        else:
            all_hits.extend(response.get('results', []))
        # Grab more hits_path if there is a next page supplied.
        next_page = response.get('nextPage', None)
        if next_page is not None and next_page < max_pages:
            next_hits = self.get_search_results(
                url,
                data,
                max_pages,
                next_page
            )
            all_hits.extend(next_hits)
        return all_hits

    def get_resp_json(self, url: str, data: Dict = None) -> Optional[Dict]:
        """
        Directly get the response json from a request.
        """
        response = requests.post(
            url,
            headers={
                'X-API-KEY': self.api_key,
                'Content-Type': 'application/json',
            },
            data=json.dumps(data),
        )
        response.raise_for_status()
        try:
            resp_json = response.json()
        except requests.exceptions.HTTPError as err:
            print(f"HTTP error: {err}")
        except requests.exceptions.ConnectionError as err:
            print(f"Connection error: {err}")
        except requests.exceptions.Timeout as err:
            print(f"Timeout error: {err}")
        except requests.exceptions.RequestException as err:
            print(f"Error: {err}")
        return resp_json

    def perform_ring_replace_search(self):
        """
        This function is used to perform a ring replacement search.
        """
        return NotImplementedError()

    def make_query(self):
        """
        This function is used to make a query to the Postera database.
        """
        pass

    def save_database_search(self):
        """
        This function is used to save the results of the Postera search as a .csv file.
        """
        pass

    def load_database_search(self):
        """
        This function is used to load the results of the Postera search from a .csv file.
        """
        pass