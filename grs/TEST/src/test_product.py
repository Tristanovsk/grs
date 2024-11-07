#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : OBS2CO
#
# ======================================================
# HISTORIQUE
# FIN-HISTORIQUE
# ======================================================


import sys
import os
import os.path
import unittest
import numpy
from grs import acutils
from grs import class_logger
from grs import product
from grs import cams_product
from grs.drivers import driver_S2_SAFE


class TestProduct(unittest.TestCase):
    """
        class for unitary test of product module
    """

    def init_log(self):
        test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = test_path + '/../output/log_file.log'
        odir = test_path + '/../output/'
        error_log_file = test_path + '/../output/error.log'
        # Init logger service (avoid to be spam by DEBUG log from rasterio...)
        class_logger.ServiceLogger(log_file=log_file, output_dir=odir, log_level='INFO', log_console=True)

    def init_prod_data(self):
        # Load data to instanciate class
        data_file = '/home/mp/penardc/scratch/PROJETS/OBS2CO/test/data/S2A_MSIL1C_20231012T104951_N0509_R051_T31TCJ_20231012T143114.SAFE'
        bandIds = range(13)
        resolution = 20
        self.l1c = driver_S2_SAFE.s2image(data_file, band_idx=bandIds, resolution=resolution)
        self.l1c.load_product()

    def test_product(self):
        """
            unitary test for product class
        """

        # Prepare log
        self.init_log()
        # Prepare test data
        self.init_prod_data()
        # instantiate product
        prod = product(self.l1c)

        self.assertEqual(prod.sensor, 'S2A')
        self.assertEqual(prod.date_str, '2023-10-12T10:49:51.024Z')
        self.assertEqual(prod.width, 5490)
        self.assertEqual(prod.height, 5490)
        self.assertAlmostEqual(prod.lonmin, 0.4959285929146376, 16)
        self.assertAlmostEqual(prod.lonmax, 1.8886830141969135, 16)
        self.assertAlmostEqual(prod.latmin, 43.238262552870076, 16)
        self.assertAlmostEqual(prod.latmax, 44.247830683507615, 16)
        self.assertAlmostEqual(prod.xmin, 300000.0, 1)
        self.assertAlmostEqual(prod.xmax, 409800.0, 1)
        self.assertAlmostEqual(prod.ymin, 4790220.0, 1)
        self.assertAlmostEqual(prod.ymax, 4900020.0, 1)
        self.assertAlmostEqual(float(prod.U), 1.0022856395562, 16)

        #load_auxiliary_data used in init, set class attributes
        # Other method not tested
        # get_flag --> used in load_flag not used
        # set_outfile --> not used
        # set_aeronetfile --> not used
        # get_elevation --> not used
        # load_flags --> not used
# class algo --> not used        
# get_elevation --> not used
