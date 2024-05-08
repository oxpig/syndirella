#!/usr/bin/env python3
"""
syndirella.cobbler.Cobbler.py

This module contains the Cobbler class. It is used to perform retrosynthetic analysis of the base compound. It handles
Poster searches, including multiple elaboration searches.
"""
import os
from typing import (List, Dict)
from syndirella.cobblers_workshop.CobblersWorkshop import CobblersWorkshop
from syndirella.Postera import Postera
from syndirella.SMARTSHandler import SMARTSHandler
from syndirella.Fairy import Fairy

class Cobbler:
    """
    A single cobbler represents 1 base compound design. It can contain multiple routes (CobblerWorkshop objects).
    """
    def __init__(self,
                 base_compound: str,
                 output_dir: str):
        self.base_compound: str = base_compound
        # Manifold API
        self.url = "https://api.postera.ai"
        self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.reaction_names = SMARTSHandler().reaction_smarts.keys()
        self.fairy = Fairy()
        self.output_dir = output_dir

    def get_routes(self) -> List[CobblersWorkshop]:
        """
        This function is used to get the routes for the base compound. Main function that is called.
        """
        # get the routes from PostEra
        routes: List[List[Dict]] = self._perform_route_search()
        # get final route(s) which will usually just be one, but will be multiple if specified in fairy filters
        final_routes: List[List[Dict]] = self._get_final_routes(routes)
        # create the cobblers workshops
        cobblers_workshops: List[CobblersWorkshop] = self._create_cobblers_workshops(final_routes)
        return cobblers_workshops

    def _perform_route_search(self, max_pages: int = 10) -> List[List[Dict]]:
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
                "vendors": ["enamine_bb", "mcule", "mcule_ultimate", 'generic']
            },
            max_pages=max_pages,
        )
        return retro_hits

    def _get_final_routes(self, routes: List[List[Dict]]) -> List[List[Dict]]:
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
            reaction_names = [reaction['name'].replace(" ","_") for reaction in reactions]
            if all([name in self.reaction_names for name in reaction_names]):
                passing_routes.append(reactions)
        try:
            # get first route and then other non-first routes if specified in fairy filters
            final_routes: List[List[Dict]] = self.fairy.get_final_routes(passing_routes)
        except IndexError:
            # there were no passing routes
            final_routes = []
            print(f"No routes found for {self.base_compound}.")
        return final_routes

    def _create_cobblers_workshops(self, final_routes: List[List[Dict]]) -> List[CobblersWorkshop]:
        """
        From the final routes, creates the cobblers workshops.
        """
        cobblers_workshops = []
        try:
            assert len(final_routes) > 0, "There are no final routes."
        except AssertionError as e:
            print(e)
            return cobblers_workshops
        for route in final_routes:
            route: List[Dict]
            # get the reactions and the product
            reaction_names: List[str] = [reaction['name'].replace(" ","_") for reaction in route]
            product = route[0]['productSmiles'] # final product
            reactants: List[str] = [reaction['reactantSmiles'] for reaction in route]
            # flip order of reaction_names and reactants to make it a forward synthesis
            reaction_names = reaction_names[::-1]
            reactants = reactants[::-1]
            self._print_route(reaction_names, reactants, product)
            # create the cobblers workshop
            cobblers_workshop = CobblersWorkshop(
                product=product,
                reactants=reactants,
                reaction_names=reaction_names,
                num_steps=len(reaction_names),
                output_dir=self.output_dir,
                filter=False,
                atoms_ids_expansion=None
            )
            cobblers_workshops.append(cobblers_workshop)
        return cobblers_workshops

    def _print_route(self, reaction_names: List[str], reactants: List[str], product: str):
        """
        This function is used to print the route.
        """
        print('Syndirella 👑 will elaborate the following route:')
        for i, reaction in enumerate(reaction_names):
            print('Step', i+1)
            print(f"{reaction} <- {reactants[i]}")
        print('Final product:', product)
        print()

    def save(self):
        """
        Pickle cobbler object, so it can be read later.
        """
        # TODO: Finish this function
        pass
