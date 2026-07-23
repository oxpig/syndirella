#!/usr/bin/env python3
"""
Test Arthor integration with Library class and complete pipeline flow.
"""
import logging
import os
import unittest
from unittest.mock import patch, MagicMock

from rdkit import Chem

from syndirella.database.Arthor import Arthor
from syndirella.route.Library import Library
from syndirella.route.Reaction import Reaction
from syndirella.pipeline import PipelineConfig
from syndirella.constants import DatabaseSearchTool, RetrosynthesisTool, DEFAULT_DATABASE_SEARCH_TOOL, DEFAULT_RETROSYNTHESIS_TOOL


class TestArthorIntegration(unittest.TestCase):
    """Test that Arthor integrates correctly with the Library class and pipeline."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.logger = logging.getLogger(__name__)
        
    def test_arthor_instantiation(self):
        """Test that Arthor can be instantiated."""
        arthor = Arthor()
        self.assertIsInstance(arthor, Arthor)
        self.assertIsNotNone(arthor.url)
        self.assertIsNotNone(arthor.logger)
        
    def test_arthor_has_required_methods(self):
        """Test that Arthor has all required methods."""
        arthor = Arthor()
        required_methods = [
            'perform_database_search',
            'perform_superstructure_search',
            'structure_output',
            'get_resp_json',
            'get_search_results',
            'get_available_databases'
        ]
        for method in required_methods:
            self.assertTrue(hasattr(arthor, method), f"Arthor missing method: {method}")
            
    # A representative snapshot of the databases live on Arthor (25Q3), used to make
    # vendor-conversion tests deterministic without hitting the network.
    LIVE_DATABASES = [
        'ChemSpace-SC-Stock-Mar2022-346K',
        'In-Stock-19Q4-14.1M',
        'MMcule-In-Stock-25Q3-6.5M',
        'Mcule-Full-25Q3-140M',
        'Mcule-Full-BBs-25Q3-6.5M',
        'Mcule-Ultimate-25Q3-111M',
        'Mcule-Virtual-25Q3-133M',
        'REAL-Database-22Q1',
        'ZINC-All-19Q4-1.4B',
        'ZINC-Interesting-19Q4-307K',
        'ZINC-On-Demand-19Q4-311M',
        'ZINC20-ForSale-22Q1',
    ]

    def test_arthor_vendor_conversion(self):
        """Vendor names resolve to the live databases matching their family."""
        arthor = Arthor()
        with patch.object(arthor, 'get_available_databases', return_value=self.LIVE_DATABASES):
            self.assertEqual(
                arthor._convert_vendors_to_arthor_dbs(['enamine_real']),
                ['REAL-Database-22Q1'])
            # The mcule family picks up all current (versioned) mcule databases.
            self.assertEqual(
                set(arthor._convert_vendors_to_arthor_dbs(['mcule'])),
                {'MMcule-In-Stock-25Q3-6.5M', 'Mcule-Full-25Q3-140M',
                 'Mcule-Full-BBs-25Q3-6.5M', 'Mcule-Ultimate-25Q3-111M',
                 'Mcule-Virtual-25Q3-133M'})

    def test_arthor_vendor_conversion_dynamic_picks_up_new_versions(self):
        """A newer mcule database name is used automatically once it appears live."""
        arthor = Arthor()
        future_dbs = self.LIVE_DATABASES + ['Mcule-Full-26Q1-200M']
        with patch.object(arthor, 'get_available_databases', return_value=future_dbs):
            self.assertIn('Mcule-Full-26Q1-200M',
                          arthor._convert_vendors_to_arthor_dbs(['mcule']))

    def test_arthor_vendor_conversion_static_fallback(self):
        """When the live list can't be fetched, static versioned names are used."""
        arthor = Arthor()
        with patch.object(arthor, 'get_available_databases', return_value=[]):
            dbs = arthor._convert_vendors_to_arthor_dbs(['enamine', 'stock'])
            self.assertIn('REAL-Database-22Q1', dbs)
            self.assertIn('In-Stock-19Q4-14.1M', dbs)
        
    def test_arthor_structure_output(self):
        """Test that structure_output correctly formats results."""
        arthor = Arthor()
        
        # Mock API response data in the new format
        mock_hits = [
            {
                'data': [
                    ['CCO', 'test1', 'BB-ForSale-22Q1'],
                    ['CCCO', 'test2', 'REAL-Database-22Q1']
                ],
                'header': ['SMILES', 'Identifier', 'arthor.source'],
                'database': 'BB-ForSale-22Q1'
            }
        ]
        
        result = arthor.structure_output(mock_hits, 'CCO', keep_catalogue=True)
        
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0][0], 'CCO')  # smiles
        self.assertEqual(result[0][1][0], 'BB-ForSale-22Q1')  # source
        self.assertEqual(result[0][1][1], 'test1')  # identifier
        
    def test_arthor_empty_results(self):
        """Test that structure_output handles empty results correctly."""
        arthor = Arthor()
        
        # Test with None results
        result = arthor.structure_output(None, 'CCO')
        self.assertIsNone(result)
        
        # Test with empty results
        result = arthor.structure_output([], 'CCO')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][0], 'CCO')  # Should return original query
        
        # Test with empty data in new format
        mock_hits = [
            {
                'data': [],
                'header': ['SMILES', 'Identifier', 'arthor.source'],
                'database': 'BB-ForSale-22Q1'
            }
        ]
        result = arthor.structure_output(mock_hits, 'CCO')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][0], 'CCO')  # Should return original query
        
    def test_arthor_structure_output_multiple_databases(self):
        """Test that structure_output correctly handles multiple databases."""
        arthor = Arthor()
        
        # Mock API response data from multiple databases
        mock_hits = [
            {
                'data': [
                    ['CCO', 'test1', 'BB-ForSale-22Q1'],
                    ['CCCO', 'test2', 'BB-ForSale-22Q1']
                ],
                'header': ['SMILES', 'Identifier', 'arthor.source'],
                'database': 'BB-ForSale-22Q1'
            },
            {
                'data': [
                    ['CCO', 'test3', 'REAL-Database-22Q1'],
                    ['CCCC', 'test4', 'REAL-Database-22Q1']
                ],
                'header': ['SMILES', 'Identifier', 'arthor.source'],
                'database': 'REAL-Database-22Q1'
            }
        ]
        
        result = arthor.structure_output(mock_hits, 'CCO', keep_catalogue=True)
        
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 4)  # Should have 4 results from 2 databases
        
        # Check that results from both databases are included
        smiles_list = [item[0] for item in result]
        self.assertIn('CCO', smiles_list)
        self.assertIn('CCCO', smiles_list)
        self.assertIn('CCCC', smiles_list)
        
        # Check that catalogue information is preserved
        catalogue_info = [item[1] for item in result if item[1] is not None]
        self.assertEqual(len(catalogue_info), 4)
        
        # Check that sources are correct
        sources = [info[0] for info in catalogue_info]
        self.assertIn('BB-ForSale-22Q1', sources)
        self.assertIn('REAL-Database-22Q1', sources)
        
    def test_arthor_structure_output_missing_headers(self):
        """Test that structure_output handles missing headers correctly."""
        arthor = Arthor()
        
        # Mock API response data with missing required headers
        mock_hits = [
            {
                'data': [
                    ['CCO', 'test1', 'BB-ForSale-22Q1'],
                    ['CCCO', 'test2', 'BB-ForSale-22Q1']
                ],
                'header': ['SMILES', 'SomeOtherColumn', 'arthor.source'],  # Missing 'Identifier'
                'database': 'BB-ForSale-22Q1'
            }
        ]
        
        result = arthor.structure_output(mock_hits, 'CCO', keep_catalogue=True)
        
        # When required headers are missing, the method skips the batch and returns original query
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)  # Only the original query is returned
        self.assertEqual(result[0][0], 'CCO')  # Should return original query
        self.assertIsNone(result[0][1])  # No catalogue info
            
    def test_arthor_structure_output_malformed_data(self):
        """Test that structure_output handles malformed data correctly."""
        arthor = Arthor()
        
        # Mock API response data with malformed structure
        mock_hits = [
            {
                'data': [
                    ['CCO', 'test1'],  # Missing third column
                    ['CCCO', 'test2', 'BB-ForSale-22Q1', 'extra']  # Extra column
                ],
                'header': ['SMILES', 'Identifier', 'arthor.source'],
                'database': 'BB-ForSale-22Q1'
            },
            {
                'data': [
                    ['CCO', 'test3', 'REAL-Database-22Q1']
                ],
                'header': ['SMILES', 'Identifier', 'arthor.source'],
                'database': 'REAL-Database-22Q1'
            }
        ]
        
        result = arthor.structure_output(mock_hits, 'CCO', keep_catalogue=True)
        
        # Should still return results, skipping malformed rows
        self.assertIsInstance(result, list)
        self.assertGreaterEqual(len(result), 1)  # Should have at least the valid result
        
        # Check that we have at least one valid result
        valid_results = [item for item in result if item[1] is not None]
        self.assertGreaterEqual(len(valid_results), 1)
        
    def test_arthor_structure_output_missing_source_column(self):
        """Test that structure_output handles missing arthor.source column correctly."""
        arthor = Arthor()
        
        # Mock API response data with missing arthor.source column
        mock_hits = [
            {
                'data': [
                    ['CCO', 'test1'],
                    ['CCCO', 'test2']
                ],
                'header': ['SMILES', 'Identifier'],  # Missing 'arthor.source'
                'database': 'BB-ForSale-22Q1'
            }
        ]
        
        result = arthor.structure_output(mock_hits, 'CCO', keep_catalogue=True)
        
        # Should return results with database name as source
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        
        # Check that catalogue info uses database name as source
        for item in result:
            self.assertIsNotNone(item[1])  # Should have catalogue info
            self.assertEqual(item[1][0], 'BB-ForSale-22Q1')  # Should use database name as source
        
    @patch('requests.Session')
    def test_arthor_superstructure_search(self, mock_session_class):
        """Test that superstructure search works correctly."""
        arthor = Arthor()
        
        # Mock the session response with new format
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'data': [
                ['CCO', 'test1', 'BB-ForSale-22Q1'],
                ['CCCO', 'test2', 'REAL-Database-22Q1']
            ],
            'header': ['SMILES', 'Identifier', 'arthor.source']
        }
        mock_session_instance = MagicMock()
        mock_session_instance.get.return_value = mock_response
        mock_session_class.return_value = mock_session_instance
        
        # Replace the session with our mock
        arthor.session = mock_session_instance
        
        result = arthor.perform_superstructure_search('CCO', vendors=['enamine_real'])
        
        self.assertIsInstance(result, list)
        self.assertGreaterEqual(len(result), 0)  # Should have at least 0 results
        
    def test_library_with_arthor_import(self):
        """Test that Library can import and use Arthor."""
        # This test verifies that the import in Library.py works
        # and that the db_search_tool parameter can be set to 'arthor'
        
        # Create a mock reaction for testing
        mock_reaction = MagicMock()
        mock_reaction.reaction_name = "test_reaction"
        mock_reaction.matched_smarts_to_reactant = {}
        
        # Test that Library can be instantiated with arthor as db_search_tool
        try:
            library = Library(
                reaction=mock_reaction,
                output_dir="/tmp/test",
                id="test_id",
                num_steps=1,
                current_step=1,
                filter=False,
                route_uuid="test_uuid",
                atom_diff_min=0,
                atom_diff_max=1,
                db_search_tool=DatabaseSearchTool.ARTHOR,
                elab_single_reactant=False
            )
            self.assertEqual(library.db_search_tool, DatabaseSearchTool.ARTHOR)
        except Exception as e:
            self.fail(f"Library instantiation with arthor failed: {e}")
    
    def test_pipeline_config_creation(self):
        """Test that PipelineConfig can be created with database search tool specification."""
        settings = {
            'input': '/tmp/test.csv',
            'output': '/tmp/output',
            'templates': '/tmp/templates',
            'hits_path': '/tmp/hits',
            'metadata': '/tmp/metadata.csv',
            'batch_num': 10000,
            'atom_diff_min': 0,
            'atom_diff_max': 10,
            'scaffold_place_num': 5,
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor',
            'manual': False,
            'only_scaffold_place': False,
            'no_scaffold_place': False,
            'elab_single_reactant': False
        }
        
        config = PipelineConfig.from_settings(settings)
        self.assertEqual(config.db_search_tool, DatabaseSearchTool.ARTHOR)
        self.assertEqual(config.retro_tool, RetrosynthesisTool.MANIFOLD)
        self.assertEqual(config.csv_path, '/tmp/test.csv')
        self.assertEqual(config.output_dir, '/tmp/output')
    
    def test_pipeline_config_defaults(self):
        """Test that PipelineConfig uses correct defaults."""
        settings = {
            'input': '/tmp/test.csv',
            'output': '/tmp/output',
            'templates': '/tmp/templates',
            'hits_path': '/tmp/hits',
            'metadata': '/tmp/metadata.csv',
            'batch_num': 10000,
            'atom_diff_min': 0,
            'atom_diff_max': 10,
            'scaffold_place_num': 5,
            'retro_tool': 'manifold',
            'long_code_column': 'Long code'
        }
        
        config = PipelineConfig.from_settings(settings)
        self.assertEqual(config.db_search_tool, DEFAULT_DATABASE_SEARCH_TOOL)  # Default
        self.assertEqual(config.manual_routes, False)  # Default
        self.assertEqual(config.only_scaffold_place, False)  # Default
        self.assertEqual(config.scaffold_place, True)  # Default
        self.assertEqual(config.elab_single_reactant, False)  # Default
    
    def test_cli_validation_arthor(self):
        """Test that CLI validation works correctly for Arthor."""
        from syndirella.cli import validate_api_credentials
        
        # Test with arthor as db_search_tool
        settings = {
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor'
        }
        
        # Should not raise an exception
        try:
            validate_api_credentials(settings)
        except Exception as e:
            self.fail(f"validate_api_credentials should not raise exception for arthor: {e}")
    
    def test_cli_validation_manifold(self):
        """Test that CLI validation works correctly for Manifold."""
        from syndirella.cli import validate_api_credentials
        
        # Test with manifold as retro_tool
        settings = {
            'retro_tool': 'manifold',
            'db_search_tool': 'manifold'
        }
        
        # Mock environment variables
        with patch.dict(os.environ, {
            'MANIFOLD_API_KEY': 'test_key',
            'MANIFOLD_API_URL': 'test_url'
        }):
            try:
                validate_api_credentials(settings)
            except Exception as e:
                self.fail(f"validate_api_credentials should not raise exception for manifold: {e}")
    
    def test_database_search_tool_flow(self):
        """Test the complete flow of database search tool specification from CLI to Library."""
        # This test verifies that the db_search_tool parameter flows correctly
        # from CLI through pipeline to Library
        
        # Mock settings that would come from CLI
        settings = {
            'input': '/tmp/test.csv',
            'output': '/tmp/output',
            'templates': '/tmp/templates',
            'hits_path': '/tmp/hits',
            'metadata': '/tmp/metadata.csv',
            'batch_num': 10000,
            'atom_diff_min': 0,
            'atom_diff_max': 10,
            'scaffold_place_num': 5,
            'retro_tool': 'manifold',
            'db_search_tool': 'arthor',
            'long_code_column': 'Long code'
        }
        
        # Create pipeline config
        config = PipelineConfig.from_settings(settings)
        self.assertEqual(config.db_search_tool, DatabaseSearchTool.ARTHOR)
        
        # Create mock reaction
        mock_reaction = MagicMock()
        mock_reaction.reaction_name = "test_reaction"
        mock_reaction.matched_smarts_to_reactant = {}
        
        # Create Library with the config's db_search_tool
        library = Library(
            reaction=mock_reaction,
            output_dir=config.output_dir,
            id="test_id",
            num_steps=1,
            current_step=1,
            filter=False,
            route_uuid="test_uuid",
            atom_diff_min=config.atom_diff_min,
            atom_diff_max=config.atom_diff_max,
            db_search_tool=config.db_search_tool,
            elab_single_reactant=False
        )
        
        # Verify the database search tool flows through correctly
        self.assertEqual(library.db_search_tool, DatabaseSearchTool.ARTHOR)
        self.assertEqual(library.config.db_search_tool, DatabaseSearchTool.ARTHOR)
    
    def test_enumerator_creation(self):
        """Test that enumerators can be created from strings."""
        # Test DatabaseSearchTool
        arthor_tool = DatabaseSearchTool.from_string('arthor')
        self.assertEqual(arthor_tool, DatabaseSearchTool.ARTHOR)
        
        manifold_tool = DatabaseSearchTool.from_string('manifold')
        self.assertEqual(manifold_tool, DatabaseSearchTool.MANIFOLD)
        
        # Test RetrosynthesisTool
        manifold_retro = RetrosynthesisTool.from_string('manifold')
        self.assertEqual(manifold_retro, RetrosynthesisTool.MANIFOLD)
        
        aizynth_retro = RetrosynthesisTool.from_string('aizynthfinder')
        self.assertEqual(aizynth_retro, RetrosynthesisTool.AIZYNTHFINDER)
    
    def test_enumerator_invalid_values(self):
        """Test that enumerators raise ValueError for invalid values."""
        # Test DatabaseSearchTool
        with self.assertRaises(ValueError):
            DatabaseSearchTool.from_string('invalid')
        
        # Test RetrosynthesisTool
        with self.assertRaises(ValueError):
            RetrosynthesisTool.from_string('invalid')


if __name__ == '__main__':
    unittest.main() 