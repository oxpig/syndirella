#!/usr/bin/env python3
"""
Arthor.py

This module contains the functionality for an Arthor search.
"""
import json
import logging
import os
import random
import sys
import time
from typing import (Any, List, Dict, Tuple, Optional)

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from rdkit import Chem

import syndirella.utils.fairy as fairy
from syndirella.database.DatabaseSearch import DatabaseSearch
from syndirella.utils.error import APIRetryLimitExceeded


class Arthor(DatabaseSearch):
    """
    This class contains information about the Arthor search. It will perform the Arthor search using the
    perform_database_search function. It will also store the results of the Arthor search as a .csv file.
    """

    def __init__(self):
        super().__init__()
        self.url = os.environ.get("ARTHOR_API_URL", "https://arthor.docking.org")
        self.api_key = os.environ.get("ARTHOR_API_KEY", None)  # Arthor may not require API key
        self.logger = logging.getLogger(f"{__name__}")
        
        # Configure session with retry strategy
        self.session = requests.Session()
        retry_strategy = requests.adapters.Retry(
            total=3,
            backoff_factor=1.0,
            status_forcelist=[500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        adapter = requests.adapters.HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

        # Cache of the databases currently available on Arthor (populated lazily).
        self._available_databases: list[str] | None = None

    # Vendor families -> the (case-insensitive) name fragments that identify their
    # databases on Arthor. Rather than pinning exact, versioned database names (which
    # change every quarter, e.g. Mcule-22Q1 -> Mcule-Full-25Q3), we discover whatever
    # databases are currently live on Arthor and match them by these fragments.
    VENDOR_NAME_FRAGMENTS: dict[str, list[str]] = {
        'mcule': ['mcule'],
        'stock': ['in-stock'],
        'zinc': ['zinc'],
        'enamine': ['enamine', 'real'],
        'real': ['real'],
    }

    # Aliases so previously accepted vendor keys keep working. They resolve to one of
    # the families above.
    VENDOR_ALIASES: dict[str, str] = {
        'enamine_real': 'enamine',
        'enamine_bb': 'enamine',
        'mcule_bb': 'mcule',
        'mcule_full': 'mcule',
        'mcule_v': 'mcule',
        'mcule_ultimate': 'mcule',
        'mcule_in_stock': 'mcule',
        'mcule_purchasable': 'mcule',
    }

    # Static fallback used only when the live Arthor database list cannot be fetched.
    # Kept current as of Arthor 25Q3 (https://arthor.docking.org).
    STATIC_VENDOR_TO_DBS: dict[str, list[str]] = {
        'mcule': ['MMcule-In-Stock-25Q3-6.5M', 'Mcule-Full-BBs-25Q3-6.5M',
                  'Mcule-Full-25Q3-140M', 'Mcule-Virtual-25Q3-133M',
                  'Mcule-Ultimate-25Q3-111M'],
        'stock': ['In-Stock-19Q4-14.1M'],
        'zinc': ['ZINC-All-19Q4-1.4B', 'ZINC-Interesting-19Q4-307K',
                 'ZINC-On-Demand-19Q4-311M', 'ZINC20-ForSale-22Q1'],
        'enamine': ['REAL-Database-22Q1'],
        'real': ['REAL-Database-22Q1'],
    }

    def perform_database_search(self,
                                reactant: Chem.Mol,
                                reaction_name: str,
                                search_type: str = "superstructure",
                                vendors: list[str] = ['enamine_real', 'stock', 'mcule', 'zinc']) -> List[str] | None:
        """
        This function is used to perform the Arthor search using the database_search_function.
        """
        if search_type != "superstructure":
            raise NotImplementedError(f"Search type '{search_type}' is not implemented. Only 'superstructure' search is supported.")
        
        # 1. Get additional similar reactant if reaction is one with additional reactants
        reactant_filters = fairy.load_reactant_filters()
        reactants: List[str] = fairy.find_similar_reactants(reactant=reactant,
                                                            reaction_name=reaction_name,
                                                            reactant_filters=reactant_filters)
        # 2. Perform the search for all
        hits_all: List[Tuple[str, Tuple[str, str] | None]] = []
        for smiles in reactants:
            if search_type == "superstructure":
                hits: List[Tuple[str, Tuple[str, str] | None]] | None = self.perform_superstructure_search(smiles, vendors=vendors)
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
        This function is used to perform the Arthor superstructure search.
        """
        self.logger.info(f"Running superstructure search for {smiles}.")
        self.logger.info(f"Querying third party services: {queryThirdPartyServices}")
        self.logger.info(f"Vendor list: {vendors}")
        if not isinstance(smiles, str):
            self.logger.error("Smiles must be a string.")
            raise TypeError("Smiles must be a string.")
        
        # Convert vendors list to Arthor database names
        arthor_dbs = self._convert_vendors_to_arthor_dbs(vendors)
        
        # Arthor uses GET requests with query parameters, not POST
        hits_all = []
        for dbname in arthor_dbs:
            try:
                self.logger.info(f"Querying Arthor database {dbname} for {smiles}.")
                response = self.session.get(
                    f'{self.url}/dt/{dbname}/search',
                    params={
                        'query': smiles,
                        'type': 'Substructure',
                        'length': 1000
                    },
                    timeout=30
                )
                
                if response.status_code == 503:
                    self.logger.warning(f'Arthor unavailable for database {dbname}')
                    continue
                    
                response.raise_for_status()
                data = response.json()
                
                if data.get('data'):
                    # data['data'] is a list of lists, and data['header'] is the list of column names
                    rows = data['data']
                    header = data.get('header')
                    if header:
                        self.logger.debug(f"Header: {header}")
                        # Store the raw data with header for processing
                        hits_all.append({
                            'data': rows,
                            'header': header,
                            'database': dbname
                        })
                    else:
                        self.logger.warning(f"No header found in response from {dbname}")

            except Exception as e:
                self.logger.error(f"Error querying database {dbname}: {str(e)}")
                continue
        
        hits_info: List[Tuple[str, Tuple[str, str] | None]] | None = self.structure_output(hits_all,
                                                                                           query_smiles=smiles,
                                                                                           keep_catalogue=keep_catalogue)

        return hits_info

    def _convert_vendors_to_arthor_dbs(self, vendors: list[str]) -> list[str]:
        """
        Convert vendor names into the Arthor database names that are actually live.

        Each vendor is resolved to a family (mcule, stock, zinc, enamine, real) and
        matched against the databases currently available on Arthor by case-insensitive
        name fragment, so the most recent versioned databases are picked up
        automatically. Passing ``'all'`` resolves to every known family. If the live
        database list cannot be fetched, a static (versioned) fallback is used.
        """
        available_databases = self.get_available_databases()
        use_dynamic = len(available_databases) > 0
        if not use_dynamic:
            self.logger.warning("Could not fetch live Arthor databases; falling back to "
                                 "static (versioned) database names.")

        # Expand 'all' to every family.
        requested_families: list[str] = []
        for vendor in vendors:
            if vendor == 'all':
                requested_families.extend(self.VENDOR_NAME_FRAGMENTS.keys())
                continue
            family = self.VENDOR_ALIASES.get(vendor, vendor)
            if family in self.VENDOR_NAME_FRAGMENTS:
                requested_families.append(family)
            else:
                self.logger.warning(f"Unknown Arthor vendor '{vendor}'; skipping.")

        arthor_dbs: list[str] = []
        for family in requested_families:
            if use_dynamic:
                fragments = self.VENDOR_NAME_FRAGMENTS[family]
                matches = [db for db in available_databases
                           if any(fragment in db.lower() for fragment in fragments)]
                if not matches:
                    self.logger.warning(f"No live Arthor databases matched vendor "
                                        f"'{family}' (fragments: {fragments}).")
            else:
                matches = self.STATIC_VENDOR_TO_DBS.get(family, [])
            for db_name in matches:
                if db_name not in arthor_dbs:
                    arthor_dbs.append(db_name)

        if not arthor_dbs:
            # Default to a commonly available in-stock database.
            fallback = [db for db in available_databases if 'in-stock' in db.lower()] \
                if use_dynamic else ['In-Stock-19Q4-14.1M']
            arthor_dbs = fallback or ['In-Stock-19Q4-14.1M']
            self.logger.warning(f"No databases resolved for vendors {vendors}; "
                                f"defaulting to {arthor_dbs}.")

        return arthor_dbs

    def get_available_databases(self, force_refresh: bool = False) -> list[str]:
        """
        Get the list of databases currently available on Arthor.

        The result is cached on the instance so repeated reactant searches do not
        re-query the ``/dt/data`` endpoint. Pass ``force_refresh=True`` to re-fetch.

        Returns
        -------
        list[str]
            List of available database names (empty if the API could not be reached).
        """
        if self._available_databases is not None and not force_refresh:
            return self._available_databases
        try:
            response = self.session.get(f'{self.url}/dt/data', timeout=30)
            response.raise_for_status()
            databases = response.json()
            self._available_databases = [db['displayName'] for db in databases]
        except Exception as e:
            self.logger.error(f"Error fetching available databases: {str(e)}")
            self._available_databases = []
        return self._available_databases
    


    def structure_output(self, 
                         hits: List[Dict] | None, 
                         query_smiles: str, 
                         keep_catalogue: bool = False) -> List[Tuple[str, Tuple[str, str] | None]] | None:
        """
        Formats output into a list of tuples with smiles and catalogue.
        """
        hits_info: List[Tuple[str, Tuple[str, str] | None]] = []
        if hits is None:
            self.logger.critical(
                f"Error with API output, returning None for superstructure search for {query_smiles}!")
            return None
        
        for hit_batch in hits:
            if not isinstance(hit_batch, dict) or 'data' not in hit_batch or 'header' not in hit_batch:
                self.logger.warning(f"Invalid hit batch structure: {hit_batch}")
                continue
                
            rows = hit_batch['data']
            header = hit_batch['header']
            database = hit_batch.get('database', 'unknown')
            
            # Find the column indices for the required fields
            try:
                smiles_idx = header.index('SMILES')
                identifier_idx = header.index('Identifier')
                source_idx = header.index('arthor.source') if 'arthor.source' in header else None
            except ValueError as e:
                self.logger.warning(f"Required column not found in header {header}: {e}")
                continue
            
            for row in rows:
                if len(row) != len(header):
                    self.logger.warning(f"Row length {len(row)} does not match header length {len(header)}: {row}")
                    continue
                    
                try:
                    smiles = row[smiles_idx]
                    identifier = row[identifier_idx] if identifier_idx < len(row) else None
                    source = row[source_idx] if source_idx is not None and source_idx < len(row) else database
                    
                    if keep_catalogue and identifier:
                        hits_info.append((smiles, (source, identifier)))
                    else:
                        hits_info.append((smiles, None))
                        
                except (IndexError, TypeError) as e:
                    self.logger.warning(f"Error processing row {row}: {e}")
                    continue
        
        if len(hits_info) == 0:
            self.logger.warning(
                f"No superstructures found for {query_smiles}, returning original query reactant, {query_smiles}.")
            # add query smiles just in case it's not returned
            hits_info.append((query_smiles, None))
        return hits_info

    @staticmethod
    def get_resp_json(url: str,
                      api_key: str = None,
                      data: Dict = None,
                      retries: int = 50,
                      backoff_factor: float = 0.5) -> Optional[Dict] | None:
        """
        Directly get the response json from a request, with retry mechanism for handling 429 status code.
        """
        logger = logging.getLogger(__name__)
        for attempt in range(retries):
            try:
                headers = {'Content-Type': 'application/json'}
                if api_key:
                    headers['X-API-KEY'] = api_key

                response = requests.post(
                    url,
                    headers=headers,
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
                        raise APIRetryLimitExceeded(smiles=data.get('query', None))

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
                           data: Dict[str, Any],
                           api_key: str = None,
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
        response: Dict | None = Arthor.get_resp_json(url, api_key, data)
        if response is None:  # API error, returning None
            return None
        elif response.get('data') is not None:
            all_hits.extend(response.get('data', []))
        elif response.get('results') is not None:
            all_hits.extend(response.get('results', []))
        # Grab more hits_path if there is a next page supplied.
        if page >= max_pages:
            return all_hits
        next_page = response.get('nextPage', None)
        if next_page is not None:
            next_hits = Arthor.get_search_results(
                url,
                data,
                api_key,
                next_page,
                max_pages
            )
            all_hits.extend(next_hits)
        return all_hits

    def make_query(self):
        """
        This function is used to make a query to the Arthor database.
        """
        pass

    def save_database_search(self):
        """
        This function is used to save the results of the Arthor search as a .csv file.
        """
        pass

    def load_database_search(self):
        """
        This function is used to load the results of the Arthor search from a .csv file.
        """
        pass


__all__ = ['Arthor']
