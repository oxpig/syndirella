#!/usr/bin/env python3
"""
syndirella.tests.test_all.py

This module contains the Slipper class. A slipper in this metaphor is the set of molecules that is the
product of a reaction.
"""

import unittest

import pandas as pd

from syndirella.cobblers_workshop._base import CobblersWorkshop
from syndirella.cobblers_workshop._cobbler_bench import CobblerBench
from ..smarts import SMARTSHandler
from ..cobblers_workshop._reaction import Reaction
from ..cobblers_workshop._library import Library
from rdkit import Chem
from syndirella.slipper.slipper_synthesizer._base import SlipperSynthesizer
from syndirella.slipper._base import Slipper
import pickle

class TestFairyFilters(unittest.TestCase):
    def setUp(self):
        self.reactants = [('COC(=O)NCCB(O)O', 'Cn1nccc1I')]
        self.product = 'COC(=O)NCCc1ccnn1C'
        self.reaction_names = ['Sp3-sp2_Suzuki_coupling']
        self.num_steps = 1
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/fairy_filters'
        self.filter = False

    def test_get_additional_reactants(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.filter)
        final_library = cobblers_workshop.get_final_library()
        self.assertIsNotNone(final_library)

class TestReactantsFiltering(unittest.TestCase):
    def setUp(self):
        self.reactants = [('OB(O)c1cccc2cc[nH]c12', 'Ic1cccc(I)n1'),
                          ('Ic1cccc(-c2cccc3cc[nH]c23)n1', 'CCCB(O)O')]
        self.product = 'CCCc1cccc(-c2cccc3cc[nH]c23)n1'
        self.reaction_names = ['Sp2-sp2_Suzuki_coupling', 'Sp3-sp2_Suzuki_coupling']
        self.num_steps = 2
        self.output_dir = \
            ('/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/'
             'reactants_filtering_w_labels')
        self.filter = False

        # need to set variables for Fragmenstein
        self.template = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/one_step_warren_A71EV2A/fragments/x0310_template.pdb'
        self.hits = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/one_step_warren_A71EV2A/fragments/clean_hits.sdf'
        self.hits_names = ['x0566_0A']
        self.n_cores = 8
        self.timeout = 240
        self.batch_num = 10000

        self.first_step_atoms_ids_expansion: (
            dict) = {0: None,
                     1: None,
                     2: True,
                     3: True,
                     4: True,
                     5: None,
                     6: None,
                     7: False,
                     8: False,
                     9: None,
                     10: None,
                     11: None,
                     12: None,
                     13: True,
                     14: None,
                     15: True}

        self.atom_ids_expansion: dict = {0: True,
                                         1: True,
                                         2: True,
                                         3: None,
                                         4: True,
                                         5: True,
                                         6: True,
                                         7: None,
                                         8: None,
                                         9: False,
                                         10: False,
                                         11: None,
                                         12: None,
                                         13: None,
                                         14: None,
                                         15: True,
                                         16: None,
                                         17: True}

    def test_get_products(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.filter)
        final_library = cobblers_workshop.get_final_library()
        slipper = Slipper(final_library,
                          atoms_ids_expansion=self.atom_ids_expansion,
                          template=self.template,
                          hits=self.hits,
                          hits_names=self.hits_names,
                          batch_num=self.batch_num)
        self.products = slipper.get_products()
        self.assertIsNotNone(self.products)
        self.placements = slipper.place_products()
        self.assertIsNotNone(self.placements)

    def test_save(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.filter)
        final_library = cobblers_workshop.get_final_library()
        final_library.save()

    def test_load(self):
        library = Library.load(self.output_dir)
        self.assertIsNotNone(library)

    def test_Delete(self):
        final_library = Library.load(self.output_dir)
        slipper = Slipper(final_library,
                          atoms_ids_expansion=self.atom_ids_expansion,
                          template=self.template,
                          hits=self.hits,
                          hits_names=self.hits_names,
                          batch_num=self.batch_num)
        slipper.output_path = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/reactants_filtering_w_labels/output'
        slipper.clean_up_placements()

    def test_GetProductsDF(self):
        final_library = Library.load(self.output_dir)
        slipper = Slipper(final_library,
                          atoms_ids_expansion=self.atom_ids_expansion,
                          template=self.template,
                          hits=self.hits,
                          hits_names=self.hits_names,
                          batch_num=self.batch_num)
        slipper.products: pd.DataFrame = (pd.read_csv('/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/reactants_filtering_w_labels/JFMKOYDGTWISRQ-UHFFFAOYSA-N_Sp3-sp2_Suzuki_coupling_products_2of2.csv'))
        slipper.output_path = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/reactants_filtering_w_labels/output'
        placements: pd.DataFrame = slipper.get_placements_df()
        self.assertIsNotNone(placements)

class TestWarrenOneStep(unittest.TestCase):
    def setUp(self):
        self.product = 'c1ccc(-c2cccc3cc[nH]c23)nc1'
        self.reactants = [('OB(O)c1cccc2cc[nH]c12', 'Brc1ccccn1')]
        self.reaction_names = ['Sp2-sp2_Suzuki_coupling']
        self.analogues = None
        self.num_steps = 1
        self.output_dir = ('/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/'
                           'one_step_warren_A71EV2A')
        self.database_search = None #postera? might not need this
        self.final_library = None
        self.slipper = None
        self.products = None

        # need to set variables for Fragmenstein
        self.template = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/one_step_warren_A71EV2A/fragments/x0310_template.pdb'
        self.hits = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/one_step_warren_A71EV2A/fragments/clean_hits.sdf'
        self.hits_names = ['x0566_0A']
        self.n_cores = 8
        self.timeout = 240
        self.batch_num = 3

        self.atom_ids_expansion: dict = {0: True,
                                         1: True,
                                         2: True,
                                         3: None,
                                         4: None,
                                         5: False,
                                         6: False,
                                         7: None,
                                         8: None,
                                         9: None,
                                         10: None,
                                         11: None,
                                         12: None,
                                         13: None,
                                         14: True}

    def test_CobblersWorkshop(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.database_search)
        self.final_library = cobblers_workshop.get_final_library()
        self.final_library.save()
        self.assertIsNotNone(self.final_library)

    def test_SlipperPlacement(self):
        self.final_library = Library.load(self.output_dir)
        self.slipper = Slipper(self.final_library, self.template, self.hits, self.hits_names, self.batch_num)
        self.products = self.slipper.get_products()
        self.placements = self.slipper.place_products()
        self.assertIsNotNone(self.placements)

    def test_AtomIDExpansionLabeling(self):
        self.final_library = Library.load(self.output_dir)
        self.slipper = Slipper(self.final_library, self.template, self.hits, self.hits_names, self.batch_num, False)
        self.products = self.slipper.get_products()
        self.assertIsNotNone(self.products)

class TestSyndirellaOneStep(unittest.TestCase):
    def setUp(self):
        self.product = 'c1cc(-c2ncco2)c2[nH]ccc2c1'
        self.reactants = [('OB(O)c1cccc2cc[nH]c12', 'Brc1ncco1')]
        self.analogues = None
        self.reaction_names = ['Sp2-sp2_Suzuki_coupling']
        self.num_steps = 1
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/syndirella/tests/test_output'
        self.database_search = None #postera? might not need this

    def test_CobblersWorkshop(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.database_search)
        final_library = cobblers_workshop.get_final_library()
        self.assertIsNotNone(final_library)

    def test_Slipper(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.database_search)
        final_library = cobblers_workshop.get_final_library()
        slipper = Slipper(final_library)
        products = slipper.get_products()
        self.assertIsNotNone(products)

class TestSyndirellaTwoStep(unittest.TestCase):
    def setUp(self):
        self.product = 'CC(C)(CCC(=O)Nc1cccc2nc(Cl)ccc12)c1nn[nH]n1'
        self.reactants = [('Nc1cccc2nc(Cl)ccc12', 'CC(C)(CCC(=O)O)B(O)O'), ('CC(C)(CCC(=O)Nc1cccc2nc(Cl)ccc12)B(O)O', 'Brc1nn[nH]n1')]
        self.reaction_names = ['Amidation', 'Sp3-sp2_Suzuki_coupling']
        self.num_steps = 2
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/syndirella/tests/test_output'
        self.database_search = None #postera? might not need this

    def test_CobblersWorkshop(self):
        cobblers_workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names, self.num_steps,
                                             self.output_dir, self.database_search)
        final_library = cobblers_workshop.get_final_library()
        self.assertIsNotNone(final_library)