# python
import logging
import os
import shutil
import sys
import unittest
from unittest.mock import patch

import syndirella.check_inputs as check_inputs
from syndirella.cli import main


class TestInputValidations(unittest.TestCase):
    def setUp(self):
        self.csv_path = '../syndirella_input_template.csv'
        self.manual_csv_path = '../syndirella_manual_input_template.csv'
        self.template_dir = 'inputs/test_inputs/templates'
        self.hits_path = 'inputs/test_inputs/A71EV2A_combined.sdf'
        self.metadata = 'inputs/test_inputs/metadata.csv'

    def test_check_csv_valid(self):
        logging.basicConfig(level=logging.INFO)
        check_inputs.check_csv(self.csv_path)
        self.assertTrue(os.path.exists(self.csv_path))

    def test_check_csv_manual_valid(self):
        logging.basicConfig(level=logging.INFO)
        check_inputs.check_csv(self.manual_csv_path)
        self.assertTrue(os.path.exists(self.manual_csv_path))

    def test_check_hit_names(self):
        logging.basicConfig(level=logging.INFO)
        check_inputs.check_hit_names(self.csv_path, self.hits_path, self.metadata, 'Long code')


class TestCLI(unittest.TestCase):
    # TODO: Need to check when getting AiZynthFinder functionality
    def setUp(self):
        self.csv_path = '../syndirella_input_template.csv'
        self.output_dir = 'outputs/test_inputs/cli'
        self.manual_csv_path = '../syndirella_manual_input_template.csv'
        self.template_dir = 'inputs/test_inputs/templates'
        self.hits_path = 'inputs/test_inputs/A71EV2A_combined.sdf'
        self.metadata = 'inputs/test_inputs/metadata.csv'
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def run_main_with_args(self, args):
        with patch.object(sys, 'argv', args):
            try:
                main()
            except SystemExit as e:
                return e.code
            return 0

    def test_valid_pipeline_arguments(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir]
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_valid_just_retro_arguments(self):
        args = ['syndirella', '--input', self.csv_path, '--output', self.output_dir, '--just_retro']
        exit_code = self.run_main_with_args(args)
        self.assertEqual(exit_code, 0)

    def test_missing_required_arguments(self):
        args = ['syndirella', '--input', self.csv_path]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                main()

    def test_missing_env_variables(self):
        args = ['syndirella', '--input', self.csv_path, '--output', 'outputs']
        with patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                main()


if __name__ == '__main__':
    unittest.main()
