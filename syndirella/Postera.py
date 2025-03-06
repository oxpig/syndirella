#!/usr/bin/env python3
"""
Postera.py

This module contains the functionality for a Postera search.
"""
import json
import logging
import os
import random
import sys
import time
from typing import (Any, List, Dict, Tuple, Optional)

import requests
from rdkit import Chem

import syndirella.fairy as fairy
from syndirella.DatabaseSearch import DatabaseSearch
from syndirella.error import APIRetryLimitExceeded


class Postera(DatabaseSearch):
    """
    This class contains information about the Postera search. It will perform the Postera search using the
    perform_database_search function. It will also store the results of the Postera search as a .csv file.
    """

    def __init__(self):
        super().__init__()
        self.url = os.environ["MANIFOLD_API_URL"]
        self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.logger = logging.getLogger(f"{__name__}")

    def perform_route_search(self,
                             compound: str) -> List[Dict[str, List[Dict[str, str]]]]:
        """
        This function is used to perform the route query.
        """
        self.logger.info(f"Running retrosynthesis analysis for {compound}...")
        if not isinstance(compound, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        retro_hits: List[Dict] | None = Postera.get_search_results(  # will be None if errored out
            url=f'{self.url}/api/v1/retrosynthesis/',
            api_key=self.api_key,
            data={
                'smiles': compound,
                "withPurchaseInfo": True,
                "vendors": ["enamine_bb", "mcule", "mcule_ultimate", 'generic']
            },
            max_pages=10
        )
        return retro_hits

    def perform_database_search(self,
                                reactant: Chem.Mol,
                                reaction_name: str,
                                search_type: str = "superstructure",
                                vendors: list[str] = ['enamine_bb', 'mcule', 'mcule_ultimate', 'enamine_real',
                                                      'enamine_made']) -> List[str] | None:
        """
        This function is used to perform the Postera search using the database_search_function.
        """
        # 1. Get additional similar reactant if reaction is one with additional reactants
        reactant_filters = fairy.load_reactant_filters()
        reactants: List[str] = fairy.find_similar_reactants(reactant=reactant,
                                                            reaction_name=reaction_name,
                                                            reactant_filters=reactant_filters)
        # 2. Perform the search for all
        hits_all: List[Tuple[str, float]] = []
        for smiles in reactants:
            if search_type == "superstructure":
                hits: List[Tuple[str, float]] | None = self.perform_superstructure_search(smiles, vendors=vendors)
                if hits is None:  # if the API query failed
                    return None
                self.logger.info(f'Found {len(hits)} hits for {smiles} before filtering.')
                hits_all.extend(hits)
        filtered_hits: List[str] = fairy.filter_molecules(hits=hits_all)
        return filtered_hits

    def perform_superstructure_search(self,
                                      smiles: str,
                                      queryThirdPartyServices: bool = False,
                                      keep_catalogue: bool = False,
                                      vendors: list[str] = ['all']) -> List[Tuple[str, Tuple[str, str] | None]] | None:
        """
        This function is used to perform the Postera superstructure search.
        """
        self.logger.info(f"Running superstructure search for {smiles}.")
        self.logger.info(f"Querying third party services: {queryThirdPartyServices}")
        self.logger.info(f"Vendor list: {vendors}")
        if not isinstance(smiles, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        superstructure_hits: List[Dict] | None = Postera.get_search_results(
            url=f'{self.url}/api/v1/superstructure/',
            api_key=self.api_key,
            data={
                'smiles': smiles,
                "entryType": "both",
                "patentDatabases": [],  # don't search over patent databases,
                "withPurchaseInfo": True,
                "queryThirdPartyServices": queryThirdPartyServices,
                "vendors": vendors
            },
            max_pages=10
        )
        hits_info: List[Tuple[str, Tuple[str, str] | None]] | None = self.structure_output(superstructure_hits,
                                                                                           query_smiles=smiles,
                                                                                           keep_catalogue=keep_catalogue)

        return hits_info

    def perform_exact_search(self,
                             smiles: str,
                             queryThirdPartyServices: bool = False,
                             catalogues: List[str] = ["all"]) -> List[Dict]:
        """
        This function is used to perform the Postera exact search.
        """
        self.logger.info(f"Running exact search for {smiles}.")
        self.logger.info(f"Searching in catalogues: {catalogues}")
        if not isinstance(smiles, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        exact_hits: List[Dict] = Postera.get_search_results(
            url=f'{self.url}/api/v1/exact/',
            api_key=self.api_key,
            data={
                'smiles': smiles,
                "queryThirdPartyServices": queryThirdPartyServices,
                "withPurchaseInfo": True,
                "vendors": catalogues
            },
            max_pages=10
        )
        return exact_hits

    def perform_exact_batch_search(self,
                                   smiles_list: List[str],
                                   queryThirdPartyServices: bool = False,
                                   catalogues: List[str] = ["all"]) -> List[Dict]:
        """
        Performs exact search with a batch of smiles.
        """
        self.logger.info(f"Running exact search for {len(smiles_list)} smiles.")
        self.logger.info(f"Searching in catalogues: {catalogues}")
        if not isinstance(smiles_list, list):
            self.logger.error("Smiles must be in a list.")
            raise TypeError("Smiles must be in a list.")
        exact_hits: List[Dict] = Postera.get_search_results(
            url=f'{self.url}/api/v1/exact/batch/',
            api_key=self.api_key,
            data={
                "smilesList": smiles_list,
                "queryThirdPartyServices": queryThirdPartyServices,
                "vendors": catalogues
            },
            max_pages=10
        )
        return exact_hits

    def structure_output(self, hits: List[Dict] | None, query_smiles: str, keep_catalogue: bool = False) -> List[
                                                                                                                Tuple[
                                                                                                                    str,
                                                                                                                    Tuple[
                                                                                                                        str, str] | None]] | None:
        """
        Formats output into a list of tuples with smiles and catalogue.
        """
        hits_info: List[Tuple[str, Tuple[str, str] | None]] = []
        if hits is None:
            self.logger.critical(
                f"Error with API output, returning empty list for superstructure search for {query_smiles}!")
            return None
        for hit in hits:
            if keep_catalogue and type(hit['catalogEntries']) is list and len(hit['catalogEntries']) > 0:
                for entry in hit['catalogEntries']:
                    entry: dict
                    hits_info.append((hit['smiles'], (entry['catalogName'], entry['catalogId'])))
            else:
                hits_info.append((hit['smiles'], None))
        if len(hits_info) == 0:
            self.logger.warning(
                f"No superstructures found for {query_smiles}, returning original query reactant, {query_smiles}.")
            # add query smiles just in case it's not returned
            hits_info.append((query_smiles, None))
        return hits_info

    @staticmethod
    def get_resp_json(url: str,
                      api_key: str,
                      data: Dict = None,
                      retries: int = 50,
                      backoff_factor: float = 0.5) -> Optional[Dict] | None:
        """
        Directly get the response json from a request, with retry mechanism for handling 429 status code.
        """
        logger = logging.getLogger(__name__)
        for attempt in range(retries):
            try:
                response = requests.post(
                    url,
                    headers={
                        'X-API-KEY': api_key,
                        'Content-Type': 'application/json',
                    },
                    data=json.dumps(data),
                )

                if response.status_code in [429, 504]:
                    if attempt < retries - 1:
                        # Calculate wait time using jittered exponential backoff strategy with at most 3 minutes
                        wait_time = backoff_factor * (2 ** attempt)
                        if wait_time > 180:
                            # choose randomly num attempts to wait for
                            wait_time = random.uniform(0, 180)
                        error_type = "Rate limit exceeded" if response.status_code == 429 else "Gateway timeout"
                        logger.warning(f"{error_type}. Waiting for {wait_time} seconds before retrying...")
                        time.sleep(wait_time)
                        continue
                    else:
                        logger.error(
                            f"Max retries exceeded with status code {response.status_code}. Please try again later. "
                            f"{response.status_code}")
                        raise APIRetryLimitExceeded(smiles=data.get('smiles', None))

                response.raise_for_status()
                return response.json()

            except requests.exceptions.HTTPError as err:
                logger.error(f"HTTP error: {err}")
            except requests.exceptions.ConnectionError as err:
                logger.error(f"Connection error: {err}")
            except requests.exceptions.Timeout as err:
                logger.error(f"Timeout error: {err}")
            except requests.exceptions.RequestException as err:
                logger.error(f"Error: {err}")
                break  # Exit the loop on non-recoverable errors

        return None

    @staticmethod
    def get_search_results(url: str,
                           api_key: str,
                           data: Dict[str, Any],
                           page: int = 1,
                           max_pages: int = 900) -> List[Dict] | None:
        """
        Recursively get all pages for the endpoint until reach null next page or
        the max_pages threshold. The default max_pages is set to 900 since the recursion limit is 1000.
        """
        current_limit = sys.getrecursionlimit()
        if max_pages > (current_limit - 100):
            raise ValueError(f"max_pages ({max_pages}) too close to recursion limit ({current_limit})")
        # List where we will gather up all the hits_path.
        all_hits = []
        data = {
            **data,
            'page': page,
        }
        response: Dict | None = Postera.get_resp_json(url, api_key, data)
        if response is None:  # API error, returning None
            return None
        elif response.get('results') is not None:
            all_hits.extend(response.get('results', []))
        elif response.get('routes') is not None:
            all_hits.extend(response.get('routes', []))
        # Grab more hits_path if there is a next page supplied.
        if page >= max_pages:
            return all_hits
        next_page = response.get('nextPage', None)
        if next_page is not None:
            next_hits = Postera.get_search_results(
                url,
                api_key,
                data,
                next_page,
                max_pages
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
