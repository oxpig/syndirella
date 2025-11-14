#!/usr/bin/env python3
"""
syndirella.route.CobblersWorkshop.py

This module contains the CobblersWorkshop class. One instance of this object is used to describe a full route.
"""

import itertools
import logging
from typing import (List, Tuple)

import hashlib
from rdkit import Chem

from syndirella.route.SMARTSHandler import SMARTSHandler
from syndirella.utils.error import SMARTSError, MolError, SingleReactantElabError
from syndirella.route.Library import Library
from syndirella.slipper.Slipper import Slipper
from .CobblerBench import CobblerBench
from syndirella.constants import DatabaseSearchTool, RetrosynthesisTool, DEFAULT_DATABASE_SEARCH_TOOL


class CobblersWorkshop():
    """
    This is the CobblersWorkshop class. It represents a full route.
    """

    def __init__(self,
                 scaffold: str,
                 reactants: List[Tuple[str]],
                 reaction_names: List[str],
                 num_steps: int,
                 output_dir: str,
                 filter: bool,
                 id: str,
                 atom_diff_min: int,
                 atom_diff_max: int,
                 elab_single_reactant: bool = False,
                 atoms_ids_expansion: dict = None,
                 db_search_tool: DatabaseSearchTool = DEFAULT_DATABASE_SEARCH_TOOL,
                 retro_tool: RetrosynthesisTool = None,
                 reference_db: str = None):
        self.logger = logging.getLogger(f"{__name__}")
        s = str((scaffold, reactants, reaction_names, num_steps))
        self.route_uuid: str = hashlib.md5(s.encode('utf-8')).hexdigest()
        self.scaffold: str = self.check_product(scaffold)
        self.id: str = id
        self.reaction_names: List[str] = self.check_reaction_names(reaction_names)
        self.reactants: List[Tuple[str]] = self.check_reactants(reactants)
        self.atom_diff_min: int = atom_diff_min
        self.atom_diff_max: int = atom_diff_max
        self.num_steps: int = num_steps
        self.elab_single_reactant_int: int | None = None
        self.elab_single_reactant: bool = self.check_and_assign_elab_setting(elab_single_reactant)
        self.output_dir: str = output_dir
        self.filter: bool = filter
        self.atoms_ids_expansion: dict = atoms_ids_expansion  # should only be internal step
        self.cobbler_benches: List[CobblerBench] = []  # To store instances for each step
        self.first_library: Library = None
        self.final_library: Library = None

        self.db_search_tool: DatabaseSearchTool = db_search_tool
        self.reference_db: str | None = reference_db
        self.retro_tool: RetrosynthesisTool = retro_tool

    def check_product(self, product: str) -> str:
        """
        Checks scaffold can be converted to a molecule and can be sanitized. If not, logs an error.
        """
        if Chem.MolFromSmiles(product) is None:
            self.logger.error(f"Could not create a molecule from the smiles {product}.")
            raise MolError(smiles=product,
                           inchi=self.id,
                           route_uuid=self.route_uuid,
                           message=f"Could not create a molecule from the smiles {product}.")
        return product

    def check_reactants(self, reactants: List[Tuple[str]]) -> List[Tuple[str]]:
        """
        Checks reactants can be converted to molecules and can be sanitized. If not, logs an error.
        """
        if self.reaction_names is None:
            self.logger.error("Reaction names not defined.")
            raise SMARTSError(message="Reaction names not defined.",
                              smiles=self.scaffold,
                              inchi=self.id,
                              route_uuid=self.route_uuid)

        for i, (reactants_step, reaction_name) in enumerate(zip(reactants, self.reaction_names)):
            if "deprotect" in reaction_name:  # Filter out None or empty reactants for deprotection steps
                reactants[i] = tuple(
                    reactant for reactant in reactants_step if reactant and reactant not in {"None", ""})

            for reactant in reactants[i]:
                if Chem.MolFromSmiles(reactant) is None:
                    self.logger.error(
                        f"Could not create a molecule from the smiles {reactant} for {reaction_name} step.")
                    raise MolError(smiles=reactant,
                                   inchi=self.id,
                                   route_uuid=self.route_uuid,
                                   message=f"Could not create a molecule from the smiles {reactant}.")

        return reactants

    def check_reaction_names(self, reaction_names: List[str]) -> List[str]:
        """
        Checks reaction names are valid and can be found in the SMARTS handler.
        """
        smarts_handler = SMARTSHandler()
        for reaction_name in reaction_names:
            if reaction_name not in smarts_handler.reaction_smarts:
                self.logger.error(f"Reaction name {reaction_name} not found in SMARTS handler.")
                raise SMARTSError(message=f"Reaction name {reaction_name} not found in SMARTS handler.",
                                  smiles=self.scaffold,
                                  inchi=self.id,
                                  route_uuid=self.route_uuid)
        return reaction_names

    def check_and_assign_elab_setting(self, elab_single_reactant: bool) -> bool:
        """
        Checks and assigns the elab_single_reactant setting. If True, assigns the reactant int to elaborate.
        """
        if elab_single_reactant:
            if self.num_steps > 1:
                self.logger.error("Single reactant elaboration only supported for single step reactions.")
                raise SingleReactantElabError(message="Single reactant elaboration only supported for single step reactions.",
                                              smiles=self.scaffold,
                                              inchi=self.id,
                                              route_uuid=self.route_uuid)
            self.elab_single_reactant_int = self.assign_single_reactant_int()
        return elab_single_reactant

    def get_cobbler_bench(self, step: int) -> CobblerBench:
        """
        This function is used to get the cobbler bench for a specific step.
        """
        reaction_name: str = self.reaction_names[step]
        reactant_smiles: Tuple[str] = self.reactants[step]
        cobbler_bench = CobblerBench(
            scaffold=self.scaffold,
            reaction_name=reaction_name,
                                      reactant_smiles=reactant_smiles,
                                      output_dir=self.output_dir,
                                      filter=self.filter,
                                      id=self.id,
                                      num_steps=self.num_steps,
                                      current_step=step+1, # 1-based
                                      route_uuid=self.route_uuid,
                                      atom_diff_min=self.atom_diff_min,
                                      atom_diff_max=self.atom_diff_max,
                                      elab_single_reactant=self.elab_single_reactant,
                                      elab_single_reactant_int=self.elab_single_reactant_int,
                                      db_search_tool=self.db_search_tool,
                                      reference_db=self.reference_db)
        return cobbler_bench

    def define_route(self):
        """
        This function is used to define the reaction and finding ids and attachment atoms.
        """
        for step in range(self.num_steps):
            cobbler_bench = self.get_cobbler_bench(step)
            self.cobbler_benches.append(cobbler_bench)

    def format_workshops_from_routes(self, routes: List[Tuple[Tuple[str, Tuple[str]], ...]]):
        """
        This function is used to format workshops from the combinations of reaction names and reactants.
        Combinations is structured like:
        [(('reaction1_name', ('reactant1_smiles','reactant2_smiles')), ('reaction2_name', ('reactant2_smiles',)), ...), # first route
        (('reaction1_name', ('reactant1_smiles','reactant2_smiles')), ('reaction2_name', ('reactant2_smiles',)), ...),...]
        """
        # can either be single route or more than one
        workshops = []
        for i, route in enumerate(routes):
            product: str = self.scaffold
            reactants: List[Tuple[str]] = [rxn[1] for rxn in route]
            reaction_names: List[str] = [rxn[0].replace(' ', '_') for rxn in route]
            num_steps: int = len(reaction_names)
            cobbler_workshop = CobblersWorkshop(scaffold=product,
                                                reactants=reactants,
                                                reaction_names=reaction_names,
                                                num_steps=num_steps,
                                                output_dir=self.output_dir,
                                                filter=self.filter,
                                                id=self.id,
                                                atoms_ids_expansion=self.atoms_ids_expansion,
                                                atom_diff_min=self.atom_diff_min,
                                                atom_diff_max=self.atom_diff_max,
                                                elab_single_reactant=self.elab_single_reactant,
                                                db_search_tool=self.db_search_tool)
            if self.elab_single_reactant:
                if reaction_names == self.reaction_names and all(
                        set(r) in [set(x) for x in self.reactants] for r in reactants):
                    # it is the original route
                    elab_int = 1 if self.elab_single_reactant_int == 0 else 0
                else:
                    elab_int = 0 if i % 2 == 0 else 1
                cobbler_workshop.elab_single_reactant_int = elab_int
            workshops.append(cobbler_workshop)
        if len(workshops) == 0:
            self.logger.error(f"No additional routes for {self.scaffold} could be created...")
        return workshops

    def configure_additional_routes(self, reaction_names_to_replace: List[str]) -> List:
        """
        This function is used to configure the additional routes from the combination of all the cobbler benches.
        Could be improved to read better...
        """
        master_benches = []
        for bench in self.cobbler_benches:
            all_benches = [rxn for rxn in bench.additional_reactions]
            master_benches.extend(all_benches)
        combinations = list(itertools.product(master_benches))
        if not self.elab_single_reactant:
            # check there is no route that matches the original if an extra route is not needed
            for route in combinations:
                if all([rxn[0] in reaction_names_to_replace for rxn in route]):
                    combinations.remove(route)
        if self.elab_single_reactant:
            # Copy each combination and insert it directly after the original
            for i in range(len(combinations)):
                combinations.insert(2 * i + 1, combinations[2 * i])
            # add the workshop reaction first in the list
            original_route = []
            for bench in self.cobbler_benches:
                original_bench = (bench.reaction_name, bench.reactant_smiles)
                original_route.append(original_bench)
            combinations.insert(0, tuple(original_route))
        additional_routes: List[CobblersWorkshop] = self.format_workshops_from_routes(combinations)
        return additional_routes

    def assign_single_reactant_int(self) -> int:
        """
        Assigns the reactant int to elaborate for the route, always assign 0 if not already assigned.
        """
        if self.elab_single_reactant_int is None:
            elab_single_reactant_int = 0
        else:
            elab_single_reactant_int = self.elab_single_reactant_int
        self.logger.info(
            f"Assigning the reactant to elaborate as {elab_single_reactant_int} out of {self.reactants[0]}")
        return elab_single_reactant_int

    def get_additional_routes(self) -> List | None:
        """
        This function is used to get an alternative route (as a CobblersWorkshop object) if it is needed such as
        containing alternative reactions or setting only one reactant to elaborate.
        """
        self.define_route()
        reaction_names_to_replace = []
        for bench in self.cobbler_benches:
            if bench.get_additional_reactions():  # Checks if alternative reaction is possible and get it
                reaction_names_to_replace.append(bench.reaction_name)
        if len(reaction_names_to_replace) == 0:
            self.logger.info(f"No alternative reactions found for {self.scaffold}.")
            if self.elab_single_reactant is False:
                return None
        additional_routes: List[CobblersWorkshop] = self.configure_additional_routes(
            reaction_names_to_replace)
        return additional_routes

    def log_route(self):
        """
        This function is used to log the route.
        """
        route_message = f"""
        
        Syndirella 👑 will elaborate the following route for {self.scaffold} | {self.id}:
        Route UUID: {self.route_uuid}
        Reaction Names: {self.reaction_names}
        Number of Steps: {self.num_steps}
        Database Search Tool: {self.db_search_tool}
        
        """
        self.logger.info(route_message)

    #################################################

    def get_final_library(self):
        """
        This function is used to get the final library of products. It dynamically handles any number of steps.
        """
        self.log_route()
        for step in range(self.num_steps):
            self.logger.info(f"Step {step + 1} in this route using {self.reaction_names[step]}")
            if len(self.cobbler_benches) == 0 or step >= len(self.cobbler_benches):
                # if no cobbler bench or step is not in cobbler benches
                cobbler_bench = self.get_cobbler_bench(step)
                self.cobbler_benches.append(cobbler_bench)
            else:
                cobbler_bench = self.cobbler_benches[step]
            current_library = cobbler_bench.find_analogues()
            if step + 1 < self.num_steps:  # you're not at the last step
                slipper = Slipper(library=current_library, atoms_ids_expansion=self.atoms_ids_expansion)
                slipper.get_products()
            # Update the final library at each step
            self.final_library = current_library
        return self.final_library
