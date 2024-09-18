#!/usr/bin/env python3
"""
syndirella.route.CobblersWorkshop.py

This module contains the CobblersWorkshop class. One instance of this object is used to describe a full route.
"""

from .CobblerBench import CobblerBench
from typing import (List, Tuple)
from syndirella.route.Library import Library
from syndirella.SMARTSHandler import SMARTSHandler
from rdkit import Chem
from syndirella.slipper.Slipper import Slipper
from syndirella.error import SMARTSError, MolError
import logging
import itertools
import shortuuid


class CobblersWorkshop():
    """
    This is the CobblersWorkshop class. It represents a full route.
    """

    def __init__(self,
                 product: str,
                 reactants: List[Tuple[str]],
                 reaction_names: List[str],
                 num_steps: int,
                 output_dir: str,
                 filter: bool,
                 id: str,
                 atom_diff_min: int,
                 atom_diff_max: int,
                 atoms_ids_expansion: dict = None):
        self.logger = logging.getLogger(f"{__name__}")
        self.route_uuid: str = shortuuid.ShortUUID().random(length=6)
        self.product: str = self.check_product(product)
        self.id: str = id
        self.reactants: List[Tuple[str]] = self.check_reactants(reactants)
        self.reaction_names: List[str] = self.check_reaction_names(reaction_names)
        self.atom_diff_min: int = atom_diff_min
        self.atom_diff_max: int = atom_diff_max

        self.num_steps: int = num_steps
        self.output_dir: str = output_dir
        self.filter: bool = filter
        self.atoms_ids_expansion: dict = atoms_ids_expansion  # should only be internal step
        self.cobbler_benches: List[CobblerBench] = []  # To store instances for each step
        self.first_library: Library = None
        self.final_library: Library = None

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
        for reactants_step in reactants:
            for reactant in reactants_step:
                if Chem.MolFromSmiles(reactant) is None:
                    self.logger.error(f"Could not create a molecule from the smiles {reactant}.")
                    raise MolError(smiles=reactant,
                                   inchi=self.id,
                                   route_uuid=self.route_uuid,
                                   message=f"Could not create a molecule from the smiles {reactant}.")
        return reactants

    def check_reaction_names(self, reaction_names: List[str]) -> List[str]:
        """
        Checks reaction names are valid. If not, logs an error.
        """
        # replace spaces with underscores if they exist
        reaction_names = [name.replace(" ", "_") for name in reaction_names]
        smarts_handler = SMARTSHandler()
        for reaction_name in reaction_names:
            if reaction_name not in smarts_handler.reaction_smarts.keys():
                self.logger.error(f"Reaction name {reaction_name} not found in SMARTS handler.")
                raise SMARTSError(message=f"Reaction name {reaction_name} not found in SMARTS handler.",
                                  smiles=self.product,
                                  inchi=self.id,
                                  route_uuid=self.route_uuid)
        return reaction_names

    def get_cobbler_bench(self, step: int) -> CobblerBench:
        """
        This function is used to get the cobbler bench for a specific step.
        """
        reactants = self.reactants[step]
        reaction_name = self.reaction_names[step]
        cobbler_bench = CobblerBench(product=self.product,
                                     reactants=reactants,
                                     reaction_name=reaction_name,
                                     output_dir=self.output_dir,
                                     id=self.id,
                                     num_steps=self.num_steps,
                                     current_step= step + 1,
                                     filter=self.filter,
                                     route_uuid=self.route_uuid,
                                     atom_diff_min=self.atom_diff_min,
                                     atom_diff_max=self.atom_diff_max)
        return cobbler_bench



    def define_reaction(self):
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
        [(('reaction1_name', ('reactant1_smiles','reactant2_smiles)), ('reaction2_name', ('reactant2_smiles',)), ...), # first route
        (('reaction1_name', ('reactant1_smiles','reactant2_smiles)), ('reaction2_name', ('reactant2_smiles',)), ...),...]
        """
        workshops = []
        for route in routes:
            product: str = self.product
            reactants: List[Tuple[str]] = [rxn[1] for rxn in route]
            reaction_names: List[str] = [rxn[0] for rxn in route]
            num_steps: int = len(reaction_names)
            cobbler_workshop = CobblersWorkshop(product=product,
                                                reactants=reactants,
                                                reaction_names=reaction_names,
                                                num_steps=num_steps,
                                                output_dir=self.output_dir,
                                                filter=self.filter,
                                                id=self.id,
                                                atoms_ids_expansion=self.atoms_ids_expansion,
                                                atom_diff_min=self.atom_diff_min,
                                                atom_diff_max=self.atom_diff_max)
            workshops.append(cobbler_workshop)
        if len(workshops) == 0:
            self.logger.error(f"No additional routes for {self.product} could be created...")
        return workshops

    def configure_additional_routes(self, reaction_names_to_replace: List[str]) -> List:
        """
        This function is used to configure the additional routes from the combination of all the cobbler benches.
        """
        master_benches = []
        for bench in self.cobbler_benches:
            if len(bench.additional_reactions) == 0:  # only add original if no additional reactions
                original_bench = (bench.reaction_name, bench.reactant_smiles)
                all_benches = [original_bench]
            else:
                all_benches = [rxn for rxn in bench.additional_reactions]
            master_benches.append(all_benches)
        combinations = list(itertools.product(*master_benches))
        # check there is no route that matches the original
        for route in combinations:
            if all([rxn[0] in reaction_names_to_replace for rxn in route]):
                combinations.remove(route)
        additional_routes: List[CobblersWorkshop] = self.format_workshops_from_routes(combinations)
        return additional_routes

    def get_additional_routes(self, edit_route: bool = True) -> List | None:
        """
        This function is used to get an alternative route (as a CobblersWorkshop object) if it is needed.

        Args:
        edit_route: bool
            If True, the reactants in the route will be directly edited via SMARTS.
            If False, an additional route containing the additional reactions specified from a Postera search will be
                returned.
        """
        self.define_reaction()
        reaction_names_to_replace = []
        for bench in self.cobbler_benches:
            if bench.get_additional_reactions():  # Checks if alternative reaction is possible and get it
                reaction_names_to_replace.append(bench.reaction_name)
        if len(reaction_names_to_replace) == 0:
            self.logger.info(f"No alternative reactions found for {self.product}.")
            return None
        additional_routes: List[CobblersWorkshop] = self.configure_additional_routes(reaction_names_to_replace)
        return additional_routes

    def log_route(self):
        """
        This function is used to log the route.
        """
        route_message = f"""
        
        Syndirella 👑 will elaborate the following route for {self.product} | {self.id}:
        Route UUID: {self.route_uuid}
        Reaction Names: {self.reaction_names}
        Number of Steps: {self.num_steps}
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



