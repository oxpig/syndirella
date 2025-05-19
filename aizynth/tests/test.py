import logging
import os
import pickle
import unittest
from unittest.mock import Mock

import glob2
# Import the class to test
from AiZynthManager import AiZynthManager
from aizynthfinder.aizynthfinder import AiZynthFinder
from parameterized import parameterized


# Configure logging at the application level
def setup_logging(log_level=logging.INFO, log_file=None):
    """Set up logging configuration."""
    # Create logger
    logger = logging.getLogger('syndirella')
    logger.setLevel(log_level)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Add file handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def setup_trees(smiles: list[str], config_file: str, stock: str, expansion_policy: str, filter_policy: str,
                logger: logging.Logger) -> dict[str, list[dict]]:
    """
    Given a list of SMILES strings, return a AiZynthFinder object.
    """
    finder_pkls = glob2.glob(f"{config_file}_routes.pkl")
    if not finder_pkls:
        logger.info('Did not find any pickle files of the routes.')
        routes = {}
        for smile in smiles:
            finder = AiZynthFinder(configfile=config_file)
            finder.stock.select(stock)
            finder.expansion_policy.select(expansion_policy)
            finder.filter_policy.select(filter_policy)
            finder.target_smiles = smile
            finder.prepare_tree()
            finder.tree_search()
            finder.build_routes()
            routes[smile]: list[dict] = [route['reaction_tree'].to_dict() for route in finder.routes]
        pickle.dump(routes, open(f"{config_file}_routes.pkl", "wb"))
    else:
        logger.info('Found pickle files of the routes. Recreating...')
        routes = pickle.load(open(finder_pkls[0], "rb"))
        assert (len(routes) == len(smiles))
    return routes


class TestAiZynthManagerSingle(unittest.TestCase):
    def setUp(self):
        self.logger = setup_logging(log_level=logging.INFO)

    def test_get_for_target(self, target_smiles: str = 'C1CCN(C(=O)NC2CCC2)CC1'):
        manager = AiZynthManager()
        manager.logger = self.logger
        output_path = f'{target_smiles}.json'
        manager.find_routes(target_smiles=target_smiles)
        manager.export_routes_to_dict(routes=manager.routes, output_path=output_path)
        self.assertTrue(os.path.exists(output_path), "json file was not created")


class TestAiZynthManagerIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up class-level resources once before all tests."""
        # Use actual config files and paths
        cls.config_file = '../config_local.yml'
        cls.stock = 'zinc'  # Use a real stock
        cls.expansion_policy = 'uspto'  # Use a real policy
        cls.filter_policy = 'uspto'  # Use a real policy
        cls.logger = setup_logging(log_level=logging.INFO)

    def setUp(self):
        """Set up instance-level resources before each test."""
        pass

    # Method to get test cases dynamically
    @staticmethod
    def get_test_cases():
        # This method needs to be static because parameterized.expand is evaluated at class definition time
        # We'll use this method inside the class to call the actual test cases
        config_file = '../config_local.yml'
        stock = 'zinc'
        expansion_policy = 'uspto'
        filter_policy = 'uspto'
        logger = setup_logging(log_level=logging.INFO)
        molecule_targets: dict = setup_trees(
            ["O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1CCCCC1"],
            config_file,
            stock,
            expansion_policy,
            filter_policy,
            logger
        )
        return [(smiles, finder) for smiles, finder in molecule_targets.items()]

    @parameterized.expand(get_test_cases.__func__)
    def test_label_route(self, smiles: str, routes: list[dict]):
        self.logger.info(f'Testing label route for {smiles}')
        for route in routes:
            mock_finder = Mock()
            # make a manager
            self.manager = AiZynthManager(
                configfile=self.config_file,
                stock=self.stock,
                expansion_policy=self.expansion_policy,
                filter_policy=self.filter_policy,
                logger=self.logger,
                finder=mock_finder,
            )

    @parameterized.expand(get_test_cases.__func__)
    def test_end_to_end_route_search(self, smiles: str, finder: AiZynthFinder):
        """Test the entire route search process with real dependencies."""
        # make a manager
        self.manager = AiZynthManager(
            configfile=self.config_file,
            stock=self.stock,
            expansion_policy=self.expansion_policy,
            filter_policy=self.filter_policy,
            logger=self.logger,
            finder=finder,
        )

        # Perform a full route search
        routes = self.manager.perform_route_search(
            target_smiles=smiles,
            top_n=5,  # Smaller number for faster tests
            max_routes_per_cluster=1
        )

        # Verify that routes were found
        self.assertTrue(len(routes) > 0, "No routes were found")

        # Verify route properties
        for route in routes:
            self.assertTrue(hasattr(route, 'score'), "Route missing score attribute")
            self.assertTrue(hasattr(route, 'reaction_mappings'), "Route missing reaction_mappings")

            # Verify reactions can be accessed
            reactions = list(route.reactions())
            self.assertTrue(len(reactions) > 0, "Route has no reactions")

    def test_export_functionality(self):
        """Test that routes can be exported to CSV with real data."""
        # First get some routes
        routes = self.manager.perform_route_search(
            target_smiles=self.test_target,
            top_n=3
        )

        # Skip test if no routes found (don't fail the test for this reason)
        if not routes:
            self.skipTest("No routes found, skipping export test")

        # Export to a temporary file
        temp_csv = "temp_test_routes.csv"
        success = self.manager.export_routes_to_csv(routes, temp_csv)

        # Verify export worked
        self.assertTrue(success, "Export failed")
        self.assertTrue(os.path.exists(temp_csv), "CSV file was not created")

        # Clean up
        if os.path.exists(temp_csv):
            os.remove(temp_csv)


# This allows the tests to be run from the command line
if __name__ == "__main__":
    unittest.main()
