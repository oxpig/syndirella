# python
import logging
import os
import shutil
import unittest

from syndirella.pipeline import run_pipeline


def handle_file_path(user_path: str) -> str:
    if os.path.isabs(user_path):
        return user_path
    return os.path.abspath(os.path.join(os.getcwd(), user_path))


class TestPipelineIntegration(unittest.TestCase):
    def setUp(self):
        self.settings = {
            'input': handle_file_path('../syndirella_input_template.csv'),
            'output': handle_file_path('outputs/test_pipeline/'),
            'templates': handle_file_path('inputs/test_inputs/templates'),
            'hits_path': handle_file_path('inputs/test_inputs/A71EV2A_combined.sdf'),
            'metadata': handle_file_path('inputs/test_inputs/metadata.csv'),
            'batch_num': 1,
            'atom_diff_min': 0,
            'atom_diff_max': 1,  # use only 1 for faster test
            'scaffold_place': True,
            'scaffold_place_num': 1,
            'long_code_column': 'Long code',
            'manual': False
        }
        if os.path.exists(self.settings['output']):
            shutil.rmtree(self.settings['output'])

    def test_pipeline_creates_output(self):
        print('This will take a while... I would get a coffee if I were you')
        logging.basicConfig(level=logging.INFO)
        run_pipeline(settings=self.settings)
        self.assertTrue(os.path.exists(self.settings['output']))
        files = os.listdir(self.settings['output'])
        self.assertGreater(len(files), 0)


if __name__ == '__main__':
    unittest.main()
