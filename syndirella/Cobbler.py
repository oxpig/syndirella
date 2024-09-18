#!/usr/bin/env python3
"""
syndirella.cobbler.Cobbler.py

This module contains the Cobbler class. It is used to perform retrosynthetic analysis of the scaffold compound. It handles
Postera searches, including multiple elaboration searches.
"""
import os
from typing import (List, Dict)
from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.Postera import Postera
from syndirella.SMARTSHandler import SMARTSHandler
from syndirella.error import NoSynthesisRoute
import logging
from rdkit import Chem
import syndirella.fairy as fairy

class Cobbler:
    """
    A single cobbler represents 1 scaffold compound design. It can contain multiple routes (CobblerWorkshop objects).
    """
    def __init__(self,
                 scaffold_compound: str,
                 output_dir: str,
                 atom_diff_min: int,
                 atom_diff_max: int):
        self.scaffold_compound: str = scaffold_compound
        self.id: str = fairy.generate_inchi_ID(self.scaffold_compound, isomeric=False)
        self.atom_diff_min: int = atom_diff_min
        self.atom_diff_max: int = atom_diff_max

        # Manifold API
        self.url = "https://api.postera.ai"
        self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.reaction_names = SMARTSHandler().reaction_smarts.keys()
        self.n_reactants_per_reaction = SMARTSHandler().n_reactants_per_reaction
        self.output_dir = output_dir
        self.logger = logging.getLogger(f"{__name__}")

    def get_routes(self) -> List[CobblersWorkshop]:
        """
        This function is used to get the routes for the scaffold compound. Main function that is called.
        """
        postera_search: Postera = Postera()
        routes: List[Dict[str, List[Dict[str, str]]]] = postera_search.perform_route_search(compound=self.scaffold_compound)
        route: CobblersWorkshop = self.get_route(routes)
        additional_routes: List[CobblersWorkshop] = route.get_additional_routes(edit_route=True)
        if additional_routes is not None:
            cobblers_workshops = [route] + additional_routes
        else:
            cobblers_workshops = [route]
        return cobblers_workshops

    def get_route(self, routes: List[Dict[str, List[Dict[str, str]]]]) -> CobblersWorkshop:
        """
        From Postera routes, choose the route, then create the cobblers workshop.
        """
        route: List[Dict] = self.choose_route(routes)
        cobblers_workshop: CobblersWorkshop = self.create_cobblers_workshop_from_Postera(route)
        return cobblers_workshop


    def create_cobblers_workshop_from_Postera(self, route: List[Dict]) -> CobblersWorkshop:
        """
        Creates a cobblers workshop from a route.
        """
        route = route[::-1] # reverse route since it is a retrosynthesis route.
        reaction_names = [reaction['name'].replace(" ", "_") for reaction in route]
        reactants = [reaction['reactantSmiles'] for reaction in route]
        product = self.scaffold_compound
        cobblers_workshop = CobblersWorkshop(
            product=product,
            reactants=reactants,
            reaction_names=reaction_names,
            num_steps=len(reaction_names),
            output_dir=self.output_dir,
            id=self.id,
            filter=False,
            atom_diff_min=self.atom_diff_min,
            atom_diff_max=self.atom_diff_max,
            atoms_ids_expansion=None
        )
        return cobblers_workshop

    def get_passing_routes(self, routes: List[Dict[str, List[Dict[str, str]]]]) -> List[List[Dict]]:
        """
        Gets the passing routes from the routes based on their reaction names are in the allowed list.
        """
        passing_routes = []
        for route in routes:
            # get the reactions
            reactions: List[Dict[str,str]] = route['reactions']
            # check if all reactions are in the reactions
            if len(reactions) == 0:
                continue
            reaction_names = [reaction['name'].replace(" ", "_") for reaction in reactions]
            if all([name in self.reaction_names for name in reaction_names]):
                passing_routes.append(reactions)
        return passing_routes

    def filter_routes(self, routes: List[List[Dict]]) -> List[List[Dict]]:
        """
        Filters routes that contain reactants with a single atom, or an intramolecular reaction that is coded as two
        reactants.
        """
        passing_routes: List[List[Dict]] = []
        for route in routes:
            reactants = [reactant for reaction in route for reactant in reaction['reactantSmiles']]
            reactant_atoms = [Chem.MolFromSmiles(reactant).GetNumAtoms() for reactant in reactants]
            if 1 in reactant_atoms: # remove routes that contain reactions with single atoms
                continue
            # reaction name and n reactants dict
            reaction_name_n_reactants = {reaction['name'].replace(" ", "_"): len(reaction['reactantSmiles']) for reaction in route}
            if any([reaction_name in self.n_reactants_per_reaction.keys() and n_reactants != self.n_reactants_per_reaction[reaction_name] for reaction_name, n_reactants in reaction_name_n_reactants.items()]):
                continue
            passing_routes.append(route)
        return passing_routes


    def choose_route(self, routes: List[Dict[str, List[Dict[str, str]]]]) -> List[Dict]:
        """
        Gets first route from List if all the reaction names are allowed and if doesn't contain a single atom as reactant.
        """
        passing_routes: List[List[Dict]] = self.get_passing_routes(routes)
        filtered_routes: List[List[Dict]] = self.filter_routes(passing_routes)
        # there were no passing routes
        final_route = filtered_routes[0] if len(filtered_routes) > 0 else []
        if len(final_route) == 0:
            self.logger.critical(f"No routes found for {self.scaffold_compound}.")
            raise NoSynthesisRoute(message=f"No routes found for {self.scaffold_compound}.",
                                   inchi=self.id,
                                   smiles=self.scaffold_compound)
        return final_route

    def _print_route(self, reaction_names: List[str], reactants: List[str], product: str):
        """
        This function is used to print the route.
        """
        self.logger.info('Syndirella 👑 will elaborate the following route:')
        for i, reaction in enumerate(reaction_names):
            self.logger.info(f'\n Step {i+1}: \n {reactants[i]} -> {reaction}')
        self.logger.info(f'Final product: {product}')

    def save(self):
        """
        Pickle cobbler object, so it can be read later.
        """
        # TODO: Finish this function
        pass
