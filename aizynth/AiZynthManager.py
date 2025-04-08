import json
import logging
import os
from collections import defaultdict
from typing import List, Dict, Optional

import numpy as np
from aizynthfinder.aizynthfinder import AiZynthFinder
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from sklearn.cluster import DBSCAN

from syndirella.SMARTSHandler import SMARTSHandler


class AiZynthManager:
    def __init__(
            self,
            config_file: str = os.getenv('AIZYNTH_CONFIG_FILE', 'config.yml'),  # default is 'config.yml'
            stock: str = 'zinc',
            expansion_policy: str = 'uspto',
            filter_policy: str = 'uspto',
            min_route_score: float = 0.5,
            clustering_eps: float = 0.2,
            clustering_min_samples: int = 2,
            logger: logging.Logger = None
    ):
        """
        Initialize the AiZynthManager.

        Args:
            config_file: Path to AiZynthFinder configuration file
            smirks_library_path: Path to SMIRKS library file
            min_route_score: Minimum score for a route to be considered
            clustering_eps: DBSCAN epsilon parameter for clustering
            clustering_min_samples: DBSCAN min_samples parameter for clustering
            logger: Logger object for logging
            finder: AiZynthFinder instance (if provided)
        """
        # Setup logging
        self.logger = logger or logging.getLogger(__name__)

        self.finder: AiZynthFinder = self._set_finder(config_file=config_file, stock=stock,
                                                      expansion_policy=expansion_policy,
                                                      filter_policy=filter_policy)

        self.min_route_score = min_route_score
        self.clustering_eps = clustering_eps
        self.clustering_min_samples = clustering_min_samples

        # Load SMIRKS library
        smirks_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'syndirella',
                                   'constants', 'RXN_SMIRKS_CONSTANTS.json')
        self.smirks_library: Dict[str, ChemicalReaction] = SMARTSHandler(rxn_smirks_path=smirks_path).reaction_smarts

    def _set_finder(self, config_file: str, stock: str, expansion_policy: str, filter_policy: str) -> AiZynthFinder:
        """
        Set the stock, expansion_policy, and filter_policy for the AiZynthManager.

        Args:
            config_file: Path to AiZynthFinder configuration file
            stock: Stock to use
            expansion_policy: Expansion policy
            filter_policy: Filter policy

        Returns:
            AiZynthFinder instance
        """
        try:
            self.logger.info(f"Loading AiZynthFinder configuration file {config_file}")
            finder = AiZynthFinder(configfile=config_file)
            finder.stock.select(stock)
            finder.expansion_policy.select(expansion_policy)
            finder.filter_policy.select(filter_policy)
        except Exception as e:
            self.logger.error("Error creating AiZynthFinder instance", e)
            raise e
        return finder

    def find_routes(self, target_smiles: str, top_n: int = 10) -> bool:
        """
        Find synthesis routes for a target molecule.

        Args:
            target_smiles: SMILES of target molecule
            top_n: Number of top routes to retrieve

        Returns:
            Boolean indicating success
        """
        try:
            self.logger.debug(f'Initializing AiZynthFinder instance')
            self.finder.target_smiles = target_smiles
            self.finder.prepare_tree()
            self.finder.tree_search()
            self.finder.build_routes()
            if not self.finder.extract_statistics()['is_solved']:
                self.logger.error(f"Could not find any routes for {target_smiles}")
                return False

            # Store top N routes
            self.routes: list[dict] = list(self.finder.routes)[:top_n]

            self.logger.info(f"Found {len(self.routes)} routes for {target_smiles}")
            return True
        except Exception as e:
            self.logger.error(f"Error finding routes: {e}")
            return False

    def filter_routes(self) -> List:
        """
        Filter routes.

        Returns:
            List of filtered routes
        """
        if not hasattr(self, 'routes'):
            self.logger.error("No routes available. Call find_routes first.")
            return []
        self.logger.info('Filtering routes...')
        self.filtered_routes = [
            route for route in self.routes
            if route['route_metadata']['is_solved'] and
               not route['reaction_tree'].is_branched()  # don't include if the route is branched
        ]

        self.logger.info(f"Filtered {len(self.routes)} routes to {len(self.filtered_routes)} routes")
        return self.filtered_routes

    def _compute_route_fingerprints(self, routes: List) -> np.ndarray:
        """
        Compute fingerprints for routes based on the reactions they contain.

        Args:
            routes: List of routes to compute fingerprints for

        Returns:
            Array of route fingerprints
        """
        # Create a set of all unique reactions across all routes
        all_reactions = set()
        for route in routes:
            for reaction in route.reactions():
                # Use reaction SMILES as identifier
                rxn_smiles = reaction.reaction_smiles()
                all_reactions.add(rxn_smiles)

        # Create a mapping from reaction to index
        reaction_to_idx = {rxn: i for i, rxn in enumerate(all_reactions)}

        # Create fingerprints for each route
        fingerprints = np.zeros((len(routes), len(all_reactions)))
        for i, route in enumerate(routes):
            for reaction in route.reactions():
                rxn_smiles = reaction.reaction_smiles()
                idx = reaction_to_idx[rxn_smiles]
                fingerprints[i, idx] = 1

        return fingerprints

    def cluster_routes(self) -> Dict[int, List]:
        """
        Cluster routes based on their reactions.

        Returns:
            Dictionary mapping cluster IDs to lists of routes
        """
        if not hasattr(self, 'filtered_routes'):
            self.logger.error("No filtered routes available. Call filter_routes_by_score first.")
            return {}

        if not self.filtered_routes:
            self.logger.warning("No routes to cluster")
            return {}

        # Compute fingerprints for routes
        fingerprints = self._compute_route_fingerprints(self.filtered_routes)

        # Cluster routes using DBSCAN
        clustering = DBSCAN(
            eps=self.clustering_eps,
            min_samples=self.clustering_min_samples,
            metric='jaccard'
        ).fit(fingerprints)

        # Group routes by cluster
        clusters = defaultdict(list)
        for i, label in enumerate(clustering.labels_):
            clusters[int(label)].append(self.filtered_routes[i])

        self.clusters = dict(clusters)

        self.logger.info(f"Clustered routes into {len(self.clusters)} clusters")
        return self.clusters

    def _match_reaction_to_smirks(self, reaction) -> Optional[str]:
        """
        Match a reaction to a SMIRKS pattern in the library.

        Args:
            reaction: Reaction object

        Returns:
            Reaction ID if a match is found, None otherwise
        """
        rxn_smiles = reaction.reaction_smiles()

        try:
            # Convert reaction SMILES to RDKit reaction object
            rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)

            # Check each SMIRKS pattern in the library
            for reaction_id, smirks in self.smirks_library.items():
                pattern = AllChem.ReactionFromSmarts(smirks)

                # Check if the reaction matches the pattern
                # This is a simplistic approach; in practice, you might need
                # more sophisticated matching logic
                if rxn.GetNumReactantTemplates() == pattern.GetNumReactantTemplates() and \
                        rxn.GetNumProductTemplates() == pattern.GetNumProductTemplates():
                    # Further check structural similarity if needed
                    # This is a placeholder for more sophisticated matching
                    return reaction_id

            return None
        except Exception as e:
            self.logger.error(f"Error matching reaction to SMIRKS: {e}")
            return None

    def assign_reactions_to_smirks(self) -> Dict[int, List]:
        """
        Assign reactions in routes to SMIRKS patterns and filter out
        routes with unassignable reactions.

        Returns:
            Dictionary mapping cluster IDs to lists of valid routes
        """
        if not hasattr(self, 'clusters'):
            self.logger.error("No clustered routes available. Call cluster_routes first.")
            return {}

        if not self.smirks_library:
            self.logger.warning("SMIRKS library is empty. All routes will be discarded.")
            return {}

        valid_clusters = defaultdict(list)

        for cluster_id, routes in self.clusters.items():
            for route in routes:
                route_valid = True

                # Store mappings from reactions to SMIRKS IDs
                route.reaction_mappings = {}

                for reaction in route.reactions():
                    reaction_id = self._match_reaction_to_smirks(reaction)

                    if reaction_id is None:
                        route_valid = False
                        break

                    route.reaction_mappings[reaction.reaction_smiles()] = reaction_id

                if route_valid:
                    valid_clusters[cluster_id].append(route)

        self.valid_clusters = dict(valid_clusters)

        # Count routes before and after filtering
        total_routes = sum(len(routes) for routes in self.clusters.values())
        valid_routes = sum(len(routes) for routes in self.valid_clusters.values())

        self.logger.info(f"Assigned reactions to SMIRKS patterns: {valid_routes}/{total_routes} routes are valid")

        return self.valid_clusters

    def get_best_routes(self, max_routes_per_cluster: int = 1) -> List:
        """
        Get the best routes from each cluster.

        Args:
            max_routes_per_cluster: Maximum number of routes to return per cluster

        Returns:
            List of best routes
        """
        if not hasattr(self, 'valid_clusters'):
            self.logger.error("No valid routes available. Call assign_reactions_to_smirks first.")
            return []

        best_routes = []

        for cluster_id, routes in self.valid_clusters.items():
            # Sort routes in cluster by score (descending)
            sorted_routes = sorted(routes, key=lambda r: r.score, reverse=True)

            # Take top N routes from each cluster
            best_routes.extend(sorted_routes[:max_routes_per_cluster])

        self.logger.info(f"Selected {len(best_routes)} best routes from {len(self.valid_clusters)} clusters")
        return best_routes

    def perform_route_search(
            self,
            target_smiles: str,
            top_n: int = 20,
            max_routes_per_cluster: int = 1
    ) -> List:
        """
        Run the full pipeline to get the best routes:
        1. Find routes
        2. Filter routes by score
        3. Cluster routes
        4. Assign reactions to SMIRKS
        5. Get best routes

        Args:
            target_smiles: SMILES of target molecule
            top_n: Number of top routes to retrieve
            max_routes_per_cluster: Maximum number of routes to return per cluster

        Returns:
            List of best routes
        """
        self.logger.info(f"Performing route search for {target_smiles} with AiZynthFinder")

        if not self.find_routes(target_smiles, top_n):
            return []

        self.filter_routes()
        self.cluster_routes()
        self.assign_reactions_to_smirks()
        best_routes = self.get_best_routes(max_routes_per_cluster)

        self.logger.info(f"Completed full pipeline for {target_smiles}. Found {len(best_routes)} best routes.")
        return best_routes

    def export_routes_to_dict(self, routes: List[dict], output_path: str | None = None) -> list | None:
        """
        Export routes to dict and can save to .json if specified.

        Args:
            routes: List of routes to export
            output_path: Path to output dict file (optional)

        Returns:
            Data to save
        """
        if output_path:
            assert '.json' in output_path, self.logger.error(
                f"Output file {output_path} does not contain .json extension")
        try:
            stored_routes = []
            for i, route in enumerate(routes):
                rxn_tree = route['reaction_tree']
                rxn_names = []
                rxn_metadata = []
                for reaction in rxn_tree.reactions():
                    rxn_names.append(reaction.smiles)
                    rxn_metadata.append(reaction.metadata)
                route_data = {
                    'smiles': self.finder.target_smiles,
                    'route_id': i,
                    'score': route['score']['state score'],
                    'num_reactions': len(list(rxn_tree.reactions())),
                    'reactions': rxn_names,
                    'metadata': rxn_metadata,
                    'smirks_mappings': '|'.join([
                        f"{r.reaction_smiles()}:{route.reaction_mappings.get(r.reaction_smiles(), 'N/A')}"
                        for r in route.reactions()
                    ]) if hasattr(route, 'reaction_mappings') else 'N/A'
                }
                stored_routes.append(route_data)
            if output_path:
                with open(output_path, 'w') as json_file:
                    json.dump(stored_routes, json_file, indent=4)
            return stored_routes
        except Exception as e:
            self.logger.error(f"Failed to export routes to dict: {e}")
            return None

    def export_reactions_to_dict(self, routes: List[dict], output_path: str | None = None) -> list | None:
        """
        Export reaction level resolution of routes to dict and can save to .json if specified.

        Args:
            routes: List of routes with reactions to export
            output_path: Path to output dict file (optional)

        Returns:
            Data to save
        """
        if output_path:
            assert '.json' in output_path, self.logger.error(
                f"Output file {output_path} does not contain .json extension")
        try:
            stored_reactions = []
            for i, route in enumerate(routes):
                rxn_tree = route['reaction_tree']
                for step, reaction in enumerate(rxn_tree.reactions()):
                    reaction_data = {
                        'smiles': self.finder.target_smiles,
                        'route_id': i,
                        'route_score': route['score']['state score'],
                        'reaction_metadata': reaction.metadata,
                        'step': step,
                        'n_steps': len(list(rxn_tree.reactions())),
                        'product': reaction.mol.smiles,
                        'reactants': tuple([umol.smiles for umol_tuple in reaction.reactants for umol in umol_tuple]),
                        'reaction': reaction.smiles,
                        'label': None
                    }
                    stored_reactions.append(reaction_data)
            if output_path:
                with open(output_path, 'w') as json_file:
                    json.dump(stored_reactions, json_file, indent=4)
            return stored_reactions
        except Exception as e:
            self.logger.error(f"Failed to export routes to dict: {e}")
            return None
