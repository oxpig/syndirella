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
import logging
import syndirella.fairy as fairy
from syndirella.DatabaseSearch import DatabaseSearch
import random


class Postera(DatabaseSearch):
    """
    This class contains information about the Postera search. It will perform the Postera search using the
    perform_database_search function. It will also store the results of the Postera search as a .csv file.
    """
    def __init__(self):
        super().__init__()
        self.url = "https://api.postera.ai"
        self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.logger = logging.getLogger(f"{__name__}")

    def perform_route_search(self,
                             compound: str,
                             max_pages: int = 10) -> List[Dict[str, List[Dict[str, str]]]]:
        """
        This function is used to perform the route query.
        """
        self.logger.info(f"Running retrosynthesis analysis for {compound}...")
        if not isinstance(compound, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        retro_hits: List[Dict] = Postera.get_search_results(
            url=f'{self.url}/api/v1/retrosynthesis/',
            api_key=self.api_key,
            data={
                'smiles': compound,
                "withPurchaseInfo": True,
                "vendors": ["enamine_bb", "mcule", "mcule_ultimate", 'generic']
            },
            max_pages=max_pages,
        )
        return retro_hits

    def perform_database_search(self,
                                reactant: Chem.Mol,
                                reaction_name: str,
                                search_type: str = "superstructure") -> Dict[str, float]:
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
                hits: List[Tuple[str, float]] = self.perform_superstructure_search(smiles)
                self.logger.info(f'Found {len(hits)} hits for {smiles} before filtering.')
                hits_all.extend(hits)
        filtered_hits: Dict[str, float] = fairy.filter_molecules(hits=hits_all)
        return filtered_hits

    def perform_superstructure_search(self,
                                      smiles: str,
                                      max_pages: int = 10) -> List[Tuple[str, float]]:
        """
        This function is used to perform the Postera superstructure search.
        """
        self.logger.info(f"Running superstructure search for {smiles}. Only searching for building blocks.")
        if not isinstance(smiles, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        superstructure_hits: List[Dict] = Postera.get_search_results(
            url=f'{self.url}/api/v1/superstructure/',
            api_key=self.api_key,
            data={
                'smiles': smiles,
                "entryType": "both",
                "patentDatabases": [], # don't search over pantent databases,
                "withPurchaseInfo": True,
                "vendors": ["all"]
            },
            max_pages=max_pages,
        )
        hits_info: List[Tuple[str, float]] = self.structure_output(superstructure_hits)
        # add query smiles just in case it's not returned
        hits_info.append((smiles, 0))
        return hits_info

    def perform_exact_search(self,
                             smiles: str,
                             queryThirdPartyServices: bool = False,
                             catalogues: List[str] = ["all"],
                             max_pages: int = 10) -> List[Dict]:
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
            max_pages=max_pages,
        )
        return exact_hits

    def perform_exact_batch_search(self,
                                   smiles_list: List[str],
                                   queryThirdPartyServices: bool = False,
                                   catalogues: List[str] = ["all"],
                                   max_pages: int = 10) -> List[Dict]:
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
            max_pages=max_pages,
        )
        return exact_hits


    def structure_output(self, hits: List[Dict]) -> List[Tuple[str, float]]:
        """
        Formats output where key is the smiles and value is the lead time. Could add other purchase info if interested.
        """
        hits_info: List[Tuple[str, float]] = []
        for hit in hits:
            # get minimum lead time possible
            # TODO: Remove lead_time from the output as it is not used.
            # lead_time = min([entry['purchaseInfo']['bbLeadTimeWeeks'] for entry in hit['catalogEntries']])
            hits_info.append((hit['smiles'], 0.0))
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
                        logger.error("Max retries exceeded. Please try again later.")
                        return None

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