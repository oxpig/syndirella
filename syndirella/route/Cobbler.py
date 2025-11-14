#!/usr/bin/env python3
"""
syndirella.cobbler.Cobbler.py

This module contains the Cobbler class. It is used to perform retrosynthetic analysis of the scaffold compound.
"""
import logging
import os
from typing import (List, Dict, Tuple)

from rdkit import Chem

import syndirella.utils.fairy as fairy

try:
    from syndirella.aizynth.AiZynthManager import AiZynthManager
except ModuleNotFoundError:
    pass

from syndirella.database.Postera import Postera
from syndirella.route.SMARTSHandler import SMARTSHandler
from syndirella.utils.error import NoSynthesisRoute, APIQueryError
from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.constants import RetrosynthesisTool, DEFAULT_RETROSYNTHESIS_TOOL, DatabaseSearchTool, DEFAULT_DATABASE_SEARCH_TOOL


class Cobbler:
    """
    A single cobbler represents 1 scaffold compound design. It can contain multiple routes (CobblerWorkshop objects).
    """

    def __init__(self,
                 scaffold_compound: str,
                 output_dir: str,
                 atom_diff_min: int,
                 atom_diff_max: int,
                 elab_single_reactant: bool,
                 retro_tool: RetrosynthesisTool = DEFAULT_RETROSYNTHESIS_TOOL,
                 db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL,
                 reference_db: str | None = None):
        self.scaffold_compound: str = scaffold_compound
        self.id: str = fairy.generate_inchi_ID(self.scaffold_compound, isomeric=False)
        self.atom_diff_min: int = atom_diff_min
        self.atom_diff_max: int = atom_diff_max
        self.elab_single_reactant: bool = elab_single_reactant

        self.retro_tool: RetrosynthesisTool = self._check_retro_tool(retro_tool)

        # Manifold API
        if self.retro_tool == RetrosynthesisTool.MANIFOLD:
            self.url = os.environ["MANIFOLD_API_URL"]
            self.api_key = os.environ["MANIFOLD_API_KEY"]
        self.reaction_names = SMARTSHandler().reaction_smarts.keys()
        self.n_reactants_per_reaction = SMARTSHandler().n_reactants_per_reaction
        self.output_dir = output_dir
        self.logger = logging.getLogger(f"{__name__}")

        self.db_search_tool: DatabaseSearchTool = db_search_tool
        self.reference_db: str | None = reference_db

    def _check_retro_tool(self, retro_tool: RetrosynthesisTool) -> RetrosynthesisTool:
        """
        Check if the specified retro tool exists.
        """
        if retro_tool not in [RetrosynthesisTool.MANIFOLD, RetrosynthesisTool.AIZYNTHFINDER]:
            raise NoSynthesisRoute(
                message=f'The specified retrosynthesis tool: {retro_tool} is not supported.',
                inchi=self.id,
                smiles=self.scaffold_compound
            )
        return retro_tool

    def _perform_retrosynthesis_search(self) -> List[Dict[str, List[Dict[str, str]]]]:
        """
        Perform the retrosynthesis search with the specified tool.
        """
        if self.retro_tool == RetrosynthesisTool.MANIFOLD:
            postera_search: Postera = Postera()
            routes: List[Dict[str, List[Dict[str, str]]]] | None = postera_search.perform_route_search(
                compound=self.scaffold_compound)
            if routes is None:  # if the API query failed
                raise APIQueryError(
                    message=f"API retrosynthesis query failed for {self.scaffold_compound}.",
                    inchi=self.id,
                    smiles=self.scaffold_compound)
        elif self.retro_tool == RetrosynthesisTool.AIZYNTHFINDER:
            # Use the default constructor since we fixed the path issue
            aizynth_search: AiZynthManager = AiZynthManager()
            
            # Create scaffold-specific directory for analysis file
            scaffold_dir = os.path.join(self.output_dir, self.id)
            os.makedirs(scaffold_dir, exist_ok=True)
            
            routes: List[Dict[str, List[Dict[str, str]]]] = aizynth_search.perform_route_search(
                target_smiles=self.scaffold_compound,
                matching_strategy='best_overall',
                save_analysis=True,
                analysis_output_path=os.path.join(scaffold_dir, f'{self.id}_aizynth_analysis.json'))
        return routes

    def get_routes(self) -> List[CobblersWorkshop]:
        """
        This function is used to get the routes for the scaffold compound. The main function that is called.
        """
        routes: List[Dict[str, List[Dict[str, str]]]] = self._perform_retrosynthesis_search()
        route: CobblersWorkshop = self.get_route(routes)
        additional_routes: List[CobblersWorkshop] = route.get_additional_routes()
        if additional_routes is not None:
            cobblers_workshops = [route] + additional_routes
        else:
            cobblers_workshops = [route]
        return cobblers_workshops

    def get_route(self, routes: List[Dict[str, List[Dict[str, str]]]]) -> CobblersWorkshop:
        """
        From routes, choose the route, then create the cobblers workshop.
        """
        route: List[Dict] = self.choose_route(routes)
        cobblers_workshop: CobblersWorkshop = self.create_cobblers_workshop_from_route(route)
        return cobblers_workshop

    def create_cobblers_workshop_from_route(self, route: List[Dict]) -> CobblersWorkshop:
        """
        This function is used to create the cobblers workshop from the route.
        """
        product: str = self.scaffold_compound
        reactants: List[Tuple[str]] = []
        reaction_names: List[str] = []
        num_steps: int = len(route)
        for i, step in enumerate(route):
            reactants.append(tuple(step['reactantSmiles']))
            reaction_names.append(step['name'].replace(' ', '_'))       
        cobblers_workshop: CobblersWorkshop = CobblersWorkshop(scaffold=product,
                                                                reactants=reactants,
                                                                reaction_names=reaction_names,
                                                                num_steps=num_steps,
                                                                output_dir=self.output_dir,
                                                                filter=False,
                                                                id=self.id,
                                                                atom_diff_min=self.atom_diff_min,
                                                                atom_diff_max=self.atom_diff_max,
                                                                elab_single_reactant=self.elab_single_reactant,
                                                                db_search_tool=self.db_search_tool,
                                                                retro_tool=self.retro_tool,
                                                                reference_db=self.reference_db)
        return cobblers_workshop

    def get_passing_routes(self, routes: List[Dict[str, List[Dict[str, str]]]]) -> List[List[Dict]]:
        """
        Gets the passing routes from the routes based on their reaction names are in the allowed list.
        """
        passing_routes = []
        for route in routes:
            # get the reactions
            reactions: List[Dict[str, str]] = route['reactions']
            # check if all reactions are in the reactions
            if len(reactions) == 0:
                continue
            reaction_names = [reaction['name'].replace(" ", "_") for reaction in reactions]
            if all([name in self.reaction_names for name in reaction_names]):
                passing_routes.append(reactions)
            else:
                not_allowed_reactions = [name for name in reaction_names if name not in self.reaction_names]
                self.logger.debug(f'Reactions in route that were not found in SMIRKS list: {not_allowed_reactions}')
        if len(passing_routes) == 0:
            raise NoSynthesisRoute(
                message=f'No routes returned contained all Syndirella reactions for {self.scaffold_compound}.',
                inchi=self.id,
                smiles=self.scaffold_compound)
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
            if 1 in reactant_atoms:  # remove routes that contain reactions with single atoms
                continue
            # reaction name and n reactants dict
            reaction_name_n_reactants = {reaction['name'].replace(" ", "_"): len(reaction['reactantSmiles']) for
                                         reaction in route}
            if any([reaction_name in self.n_reactants_per_reaction.keys() and n_reactants !=
                    self.n_reactants_per_reaction[reaction_name] for reaction_name, n_reactants in
                    reaction_name_n_reactants.items()]):
                continue
            last_product: str = [reaction['productSmiles'] for reaction in route][-1]
            if fairy.generate_inchi_ID(last_product, isomeric=False) != self.id: # last product isn't correct
                continue
            passing_routes.append(route)
        if len(passing_routes) == 0:
            raise NoSynthesisRoute(
                message=f'Filtered out all routes due to single atom reactants or intramolecular reactions for {self.scaffold_compound}.',
                inchi=self.id,
                smiles=self.scaffold_compound)
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
            self.logger.info(f'\n Step {i + 1}: \n {reactants[i]} -> {reaction}')
        self.logger.info(f'Final product: {product}')
