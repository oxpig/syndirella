#!/usr/bin/env python3
"""
syndirella.tests.hippo.hippo_tests.py
"""

import unittest

import pandas as pd
import pickle
from rdkit import Chem

from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.Cobbler import Cobbler
from syndirella.route.Library import Library
from syndirella.slipper.Slipper import Slipper
from syndirella.slipper.SlipperFitter import SlipperFitter
from syndirella.pipeline import run_pipeline

class TestBasicPipeline(unittest.TestCase):
    def setUp(self):
        self.test_csv = "/Users/kate_fieseler/PycharmProjects/EV-A71-2A-syndirella-run/batches_run3/batch3.csv"
        self.template_path = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/hippo/x0310_relaxed_apo.pdb'
        self.hits_path = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/hippo/A71EV2A/A71EV2A_combined.sdf'
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/batch1'
        self.batch_num = 10
        self.additional_info = ['compound_set']
        self.manual_routes = True

    def test_pipeline(self):
        run_pipeline(csv_path=self.test_csv, output_dir=self.output_dir, template_dir='/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/hippo/', hits_path=self.hits_path,
                     metadata_path=str, batch_num=self.batch_num, additional_columns=self.additional_info,
                     manual_routes=self.manual_routes)

class TestFromSlipper(unittest.TestCase):
    def setUp(self):
        self.library = pickle.load(open('/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/hippo/JHZMHCCPHXTVDH-UHFFFAOYSA-N/test_final_library.pkl', 'rb'))
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/hippo/JHZMHCCPHXTVDH-UHFFFAOYSA-N'
        self.placements = pd.read_pickle('/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/hippo/JHZMHCCPHXTVDH-UHFFFAOYSA-N/JHZMHCCPHXTVDH-UHFFFAOYSA-N_Williamson_ether_synthesis_products_3of3_placements.pkl.gz')
        self.uuid = '4LAMhz'
        self.batch_num = 20

    def test_from_slipper(self):
        slipper = Slipper(library=self.library)
        slipper.placements = self.placements
        slipper.batch_num = self.batch_num
        slipper.write_products_to_hippo(uuid=self.uuid)
