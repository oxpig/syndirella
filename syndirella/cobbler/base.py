#!/usr/bin/env python3
"""
syndirella.cobbler.base.py

This module contains the Cobbler class. It is used to perform retrosynthetic analysis of the base compound. It handles
Poster searches, including multiple elaboration searches.
"""
import os
from typing import (List, Dict, Tuple, Union, Optional)
from syndirella.cobblers_workshop._base import CobblersWorkshop
from syndirella.cobblers_workshop._postera import Postera
from syndirella.smarts import SMARTSHandler

class Cobbler:
    """
    A single cobbler represents 1 base compound design. It can contain multiple routes (CobblerWorkshop objects).
    """
    def __init__(self,
               base_compound: str):
        self.base_compound: str = base_compound
        # Manifold API
        self.url = "https://api.postera.ai"
        self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.reaction_names = SMARTSHandler().reaction_smarts.keys()

    def get_routes(self) -> List[CobblersWorkshop]:
        """
        This function is used to get the routes for the base compound.
        """
        # get the routes
        routes: List[CobblersWorkshop] = self._get_routes_from_postera()
        # get the final libraries for each route
        for route in routes:
            route.get_final_library()
        return routes

    def _get_routes_from_postera(self) -> List[Dict]:
        """
        This function is used to get the routes from PostEra.
        """
        # get the routes from PostEra
        routes = self._perform_route_search()
        # get final route(s) which will usually just be one, but will be multiple if specified in fairy filters
        final_routes = self._get_final_routes(routes)
        return final_routes

    def _perform_route_search(self, max_pages: int = 10) -> List[Dict]:
        """
        This function is used to perform the route query.
        """
        print(f"Running retrosynthesis analysis for {self.base_compound}...")
        assert type(self.base_compound) == str, "Smiles must be a string."
        retro_hits: List[Dict] = Postera.get_search_results(
            url=f'{self.url}/api/v1/retrosynthesis/',
            api_key=self.api_key,
            data={
                'smiles': self.base_compound,
                "withPurchaseInfo": True,
                #"vendors": ['all']
                "vendors": ["enamine_bb", "mcule", "mcule_ultimate", 'generic']
            },
            max_pages=max_pages,
        )
        return retro_hits

    def _get_final_routes(self, routes: List[Dict]) -> List[Dict]:
        """
        This function is used to get the final routes, which are routes that contain all reactions we have encoded.
        Then the final routes is either the first 1 or other specified by the fairy filters.
        """
        passing_routes = []
        for route in routes:
            # get the reactions
            reactions = route['reactions']
            # check if all reactions are in the reactions
            if len(reactions) == 0:
                continue
            reactions = [reaction['name'].replace(" ","_") for reaction in reactions]
            if all([reaction in self.reaction_names for reaction in reactions]):
                passing_routes.append(route)
        return passing_routes