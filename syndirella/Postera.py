#!/usr/bin/env python3
"""
Postera.py

This module contains the functionality for a Postera search.
"""
import os
from typing import (Any, List, Dict, Tuple, Optional)
from rdkit import Chem
import requests
import json
import time
from syndirella.Fairy import Fairy
from syndirella.DatabaseSearch import DatabaseSearch
from syndirella.cobblers_workshop.Reaction import Reaction


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
        self.fairy = Fairy()

    def perform_database_search(self, reactant: Chem.Mol):
        """
        This function is used to perform the Postera search using the database_search_function.
        """
        # 1. Get additional similar reactant if reaction is one with additional reactants
        reactants: List[str] = self.fairy.find(reactant, self.reaction.reaction_name)
        # 2. Perform the search for all
        hits_all: List[Tuple[str, float]] = []
        for smiles in reactants:
            if self.search_type == "superstructure":
                hits: List[Tuple[str, float]] = self.perform_superstructure_search(smiles)
                print(f'Found {len(hits)} for {smiles} before filtering.')
                hits_all.extend(hits)
        filtered_hits: Dict[str, float] = self.fairy.filter(hits_all)
        return filtered_hits

    def perform_superstructure_search(self,
                                      smiles: str,
                                      max_pages: int = 10) -> List[Tuple[str, float]]:
        """
        This function is used to perform the Postera superstructure search.
        """
        print(f"Running superstructure search for {smiles}. Only searching for building blocks.")
        assert type(smiles) == str, "Smiles must be a string."
        superstructure_hits: List[Dict] = Postera.get_search_results(
            url=f'{self.url}/api/v1/superstructure/',
            api_key=self.api_key,
            data={
                'smiles': smiles,
                "entryType": "building_block",
                "withPurchaseInfo": True,
                "vendors": ["enamine_bb", "mcule", "mcule_ultimate"]
            },
            max_pages=max_pages,
        )
        hits_info: List[Tuple[str, float]] = self.structure_output(superstructure_hits)
        # add query smiles just in case it's not returned
        hits_info.append((smiles, 0))
        return hits_info

    def structure_output(self, hits: List[Dict]) -> List[Tuple[str, float]]:
        """
        Formats output where key is the smiles and value is the lead time. Could add other purchase info if interested.
        """
        hits_info: List[Tuple[str, float]] = []
        for hit in hits:
            # get minimum lead time possible
            lead_time = min([entry['purchaseInfo']['bbLeadTimeWeeks'] for entry in hit['catalogEntries']])
            hits_info.append((hit['smiles'], lead_time))
            # could do price but seems like too much faff
        return hits_info

    @staticmethod
    def get_resp_json(url: str,
                      api_key: str,
                      data: Dict = None,
                      retries: int = 50,
                      backoff_factor: float = 0.5) -> Optional[Dict]:
        """
        Directly get the response json from a request, with retry mechanism for handling 429 status code.
        """
        for attempt in range(retries):
            response = requests.post(
                url,
                headers={
                    'X-API-KEY': api_key,
                    'Content-Type': 'application/json',
                },
                data=json.dumps(data),
            )
            if response.status_code == 429:
                if attempt < retries - 1:
                    # Calculate wait time using exponential backoff strategy
                    wait_time = backoff_factor * (2 ** attempt)
                    print(f"Rate limit exceeded. Waiting for {wait_time} seconds before retrying...")
                    time.sleep(wait_time)
                    continue
                else:
                    print("Max retries exceeded. Please try again later.")
                    return None
            response.raise_for_status()
            try:
                return response.json()
            except requests.exceptions.HTTPError as err:
                print(f"HTTP error: {err}")
            except requests.exceptions.ConnectionError as err:
                print(f"Connection error: {err}")
            except requests.exceptions.Timeout as err:
                print(f"Timeout error: {err}")
            except requests.exceptions.RequestException as err:
                print(f"Error: {err}")
            break  # If the request was successful, break out of the loop.

    # @staticmethod
    # def get_resp_json(url: str,
    #                   api_key: str,
    #                   data: Dict = None) -> Optional[Dict]:
    #     """
    #     Directly get the response json from a request.
    #     """
    #     response = requests.post(
    #         url,
    #         headers={
    #             'X-API-KEY': api_key,
    #             'Content-Type': 'application/json',
    #         },
    #         data=json.dumps(data),
    #     )
    #     response.raise_for_status()
    #     try:
    #         resp_json = response.json()
    #     except requests.exceptions.HTTPError as err:
    #         print(f"HTTP error: {err}")
    #     except requests.exceptions.ConnectionError as err:
    #         print(f"Connection error: {err}")
    #     except requests.exceptions.Timeout as err:
    #         print(f"Timeout error: {err}")
    #     except requests.exceptions.RequestException as err:
    #         print(f"Error: {err}")
    #     return resp_json

    @staticmethod
    def get_search_results(url: str,
                           api_key: str,
                           data: Dict[str, Any],
                           max_pages: int = 5,
                           page: int = 1) -> List[Dict]:
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
        response = Postera.get_resp_json(url, api_key, data)
        if response is None:
            return all_hits
        elif response.get('results') is not None:
            all_hits.extend(response.get('results', []))
        elif response.get('routes') is not None:
            all_hits.extend(response.get('routes', []))
        # Grab more hits_path if there is a next page supplied.
        next_page = response.get('nextPage', None)
        if next_page is not None and next_page < max_pages:
            next_hits = Postera.get_search_results(
                url,
                api_key,
                data,
                max_pages,
                next_page
            )
            all_hits.extend(next_hits)
        return all_hits

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