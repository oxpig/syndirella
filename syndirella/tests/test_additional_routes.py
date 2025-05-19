# python
import logging
import os
import shutil
import unittest

import pandas as pd

from syndirella.Cobbler import Cobbler
from syndirella.fairy import generate_inchi_ID
from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.route.Library import Library
from syndirella.slipper.Slipper import Slipper


def handle_file_path(user_path):
    if os.path.isabs(user_path):
        return user_path
    return os.path.abspath(os.path.join(os.getcwd(), user_path))


class TestAdditionalRouteGeneration(unittest.TestCase):
    def setUp(self):
        self.product = 'CC(=O)Nc1cncc(CC(=O)NC(C)=O)c1'
        self.reactants = [('CCOC(=O)Cc1cncc(N)c1', 'CC(=O)Cl')]
        self.reaction_names = ['Amide_Schotten-Baumann_with_amine']
        self.num_steps = 1
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/additional_route'
        self.filter = False
        self.id = generate_inchi_ID(self.product)
        self.elab_single_reactant = False
        self.atom_diff_min = 0
        self.atom_diff_max = 10

    def test_get_additional_routes(self):
        logging.basicConfig(level=logging.INFO)
        workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names,
                                    self.num_steps, self.output_dir, self.filter, id=self.id,
                                    elab_single_reactant=self.elab_single_reactant, atom_diff_max=self.atom_diff_max,
                                    atom_diff_min=self.atom_diff_min)
        routes = workshop.get_additional_routes(edit_route=True)
        self.assertIsInstance(routes, list)


class TestElabSingleReactant(unittest.TestCase):
    def setUp(self):
        self.product = 'CC(=O)Nc1cncc(CC(=O)NC(C)=O)c1'
        self.reactants = [('CCOC(=O)Cc1cncc(N)c1', 'CC(=O)Cl')]
        self.reaction_names = ['Amide_Schotten-Baumann_with_amine']
        self.num_steps = 1
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/additional_route'
        self.filter = False
        self.id = generate_inchi_ID(self.product)
        self.elab_single_reactant = True  ###
        self.atom_diff_min = 0
        self.atom_diff_max = 10
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        # self.template_path =
        # self.hits_path =
        # self.hits =
        # self.batch_num =
        # self.additional_info =
        # self.scaffold_placements =

    def test_get_additional_routes(self):
        logging.basicConfig(level=logging.INFO)
        workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names,
                                    self.num_steps, self.output_dir, self.filter, id=self.id,
                                    elab_single_reactant=self.elab_single_reactant, atom_diff_min=self.atom_diff_min,
                                    atom_diff_max=self.atom_diff_max)
        routes = workshop.get_additional_routes(edit_route=True)
        self.assertIsInstance(routes, list)

    def test_get_final_library(self):
        logging.basicConfig(level=logging.INFO)
        workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names,
                                    self.num_steps, self.output_dir, self.filter, id=self.id,
                                    elab_single_reactant=self.elab_single_reactant, atom_diff_min=self.atom_diff_min,
                                    atom_diff_max=self.atom_diff_max)
        final_library = workshop.get_final_library()
        slipper = Slipper(library=final_library)
        slipper.get_products()
        self.assertIsInstance(final_library, Library)

    def test_len_of_reactant_dataframes(self):
        logging.basicConfig(level=logging.INFO)
        workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names,
                                    self.num_steps, self.output_dir, self.filter, id=self.id,
                                    elab_single_reactant=self.elab_single_reactant, atom_diff_min=self.atom_diff_min,
                                    atom_diff_max=self.atom_diff_max)
        cobbler_bench = workshop.get_cobbler_bench(0)  # first and only step
        current_library = cobbler_bench.find_analogues()
        df1 = current_library.analogues_dataframes['r1'][0]
        df2 = current_library.analogues_dataframes['r2'][0]
        self.assertTrue((len(df1) == 1 and len(df2) > 1) or (len(df1) > 1 and len(df2) == 1))


class TestAdditionalRouteandSingleElab(unittest.TestCase):
    def setUp(self):
        self.product = 'CC(=O)Nc1cncc(CC(=O)NC(C)=O)c1'
        self.reactants = [('CCOC(=O)Cc1cncc(N)c1', 'CC(=O)Cl')]
        self.reaction_names = ['Amide_Schotten-Baumann_with_amine']
        self.num_steps = 1
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/additional_route'
        self.filter = False
        self.id = generate_inchi_ID(self.product)
        self.elab_single_reactant = True  ###
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        self.atom_diff_min = 0
        self.atom_diff_max = 10

    def test_get_additional_routes(self):
        logging.basicConfig(level=logging.INFO)
        workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names,
                                    self.num_steps, self.output_dir, self.filter, id=self.id,
                                    elab_single_reactant=self.elab_single_reactant, atom_diff_min=self.atom_diff_min,
                                    atom_diff_max=self.atom_diff_max)
        routes = workshop.get_additional_routes(edit_route=True)
        self.assertIsInstance(routes, list)
        self.assertTrue(len(routes) == 3)  # one route for the original and two for the new
        route0 = workshop
        route1 = routes[0]
        route2 = routes[1]
        route3 = routes[2]
        self.assertTrue((route0.reaction_names[0] == 'Amide_Schotten-Baumann_with_amine') and
                        (route1.reaction_names[0] == 'Amide_Schotten-Baumann_with_amine') and
                        (route2.reaction_names[0] == 'Amidation') and
                        (route3.reaction_names[0] == 'Amidation'))
        self.assertTrue((route0.elab_single_reactant_int != route1.elab_single_reactant_int) and
                        (route2.elab_single_reactant_int != route3.elab_single_reactant_int))


class TestAutoSingleElabHippo(unittest.TestCase):
    def setUp(self):
        self.output_csv_path = 'outputs/test_additional_routes/single_elab/output.csv'
        self.product = 'COCC(=O)N[C@H]1CCCCC12CCCCC2'
        self.output_dir = 'outputs/test_additional_routes/single_elab'
        self.hits_names = ['CHIKV_MacB-x1444_A_401_CHIKV_MacB-x0300+A+401+1',
                           'CHIKV_MacB-x0692_C_401_CHIKV_MacB-x0300+A+401+1']
        self.hits_path = 'inputs/test_error_handling/CHIKV_Mac_combined.sdf'
        self.template_path = handle_file_path('inputs/test_error_handling/cx0270a_apo-desolv.pdb')
        self.additional_info = {'compound_set': 'test'}
        self.csv_path = 'inputs/test_error_handling/single_elab.csv'
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_full_hippo(self):
        logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler()])
        cobbler = Cobbler(scaffold_compound=self.product,
                          output_dir=self.output_dir,
                          atom_diff_min=0,
                          atom_diff_max=3,
                          elab_single_reactant=True)
        cobbler_workshops = cobbler.get_routes()
        for workshop in cobbler_workshops:
            final_library = workshop.get_final_library()
            slipper = Slipper(library=final_library, template=self.template_path, hits_path=self.hits_path,
                              hits_names=self.hits_names,
                              batch_num=3, atoms_ids_expansion=None, additional_info=None,
                              scaffold_placements=None)
            slipper.get_products()
            slipper.place_products()
            hippo_path: str = slipper.write_products_to_hippo()  # only write at the end after placement, to get correct route_uuid
            slipper.clean_up_placements()
            self.assertTrue(os.path.exists(hippo_path))
            hippo_df = pd.read_pickle(hippo_path)
            self.assertTrue(hippo_df.shape[0] > 0)
            self.assertIn('1_single_reactant_elab', hippo_df.columns)
            self.assertTrue(hippo_df['1_single_reactant_elab'].all())


if __name__ == '__main__':
    unittest.main()
