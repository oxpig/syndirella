# python
import logging
import os
import shutil
import unittest

from syndirella.error import NoSynthesisRoute
from syndirella.pipeline import elaborate_compound_full_auto


def handle_file_path(user_path):
    if os.path.isabs(user_path):
        return user_path
    return os.path.abspath(os.path.join(os.getcwd(), user_path))


class TestErrorHandlingNoSynthesisRoute(unittest.TestCase):
    # TODO: Need to check when getting AiZynthFinder functionality
    def setUp(self):
        self.output_csv_path = 'outputs/test_error_handling/no_synthesis_route/output.csv'
        self.product = 'Cn1cc(C(=O)NCS(=O)(=O)c2ncc[nH]2)c2ccccc21'
        self.output_dir = 'outputs/test_error_handling/no_synthesis_route'
        self.hits_names = ['CHIKV_MacB-x0270_A_304_CHIKV_MacB-x0300+A+401+1',
                           'CHIKV_MacB-x1268_B_304_CHIKV_MacB-x0300+A+401+1']
        self.hits_path = 'inputs/test_error_handling/no_synthesis_route/CHIKV_Mac_combined.sdf'
        self.template_path = handle_file_path('inputs/test_error_handling/no_synthesis_route/cx0270a_apo-desolv.pdb')
        self.additional_info = {'compound_set': 'test'}
        self.csv_path = 'inputs/test_error_handling/no_synthesis_route/test.csv'

    def test_no_synthesis_route_error(self):
        logging.basicConfig(level=logging.INFO)
        # remove output directory if it exists
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        try:
            elaborate_compound_full_auto(product=self.product,
                                         output_dir=self.output_dir,
                                         hits=self.hits_names,
                                         hits_path=self.hits_path,
                                         batch_num=3,
                                         template_path=self.template_path,
                                         additional_info=self.additional_info,
                                         csv_path=self.csv_path,
                                         atom_diff_min=0,
                                         atom_diff_max=10,
                                         scaffold_place_num=1,
                                         only_scaffold_place=False,
                                         scaffold_place=False
                                         )
        except NoSynthesisRoute as e:
            with open(self.output_csv_path, 'w') as f:
                f.write('Error,Message\n')
                f.write(f'NoSynthesisRoute,{str(e)}\n')
        self.assertTrue(os.path.exists(self.output_csv_path))
        self.assertGreater(os.path.getsize(self.output_csv_path), 0)


if __name__ == '__main__':
    unittest.main()
