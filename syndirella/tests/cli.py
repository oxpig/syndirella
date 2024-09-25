import unittest
import subprocess
from syndirella.pipeline import run_pipeline

class TestCLI(unittest.TestCase):
    def test_cli(self):
        subprocess.check_call(['syndirella', '--help'])

    def test_crystal_name_cli(self):
        subprocess.check_call(['syndirella',
                               '--input', '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/designs/D68EV3CPROA_knitwork_722_syndirella_master.csv',
                               '--output', 'test',
                               '--templates', '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/fragments/templates',
                               '--hits_path', '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/fragments/D68EV3CPROA_combined.sdf',
                               '--metadata', '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/fragments/metadata.csv'])

    def setUp(self):
        self.settings = {
            'input': '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/designs/D68EV3CPROA_knitwork_722_syndirella_master.csv',
            'output': 'test',
            'templates': '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/fragments/templates',
            'hits_path': '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/fragments/D68EV3CPROA_combined.sdf',
            'metadata': '/Users/kate_fieseler/PycharmProjects/EV-D68-3C-PROA-syndirella-run/fragments/metadata.csv',
            'batch_num': 5,
            'atom_diff_min': 0,
            'atom_diff_max': 10,
            'scaffold_place_num': 3
        }

    def test_crystal_name(self):
        run_pipeline(self.settings)