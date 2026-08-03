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
from pathlib import Path
from grs import acutils
from grs import class_logger
from grs.product import Product
import xarray as xr


class TestProduct(unittest.TestCase):
    """
        class for unitary test of product module
    """

    test_path = ""

    @classmethod
    def setUpClass(cls) -> None:
        cls.test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = cls.test_path + '/../output/log_file.log'
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        odir = cls.test_path + '/../output/'
        logger = class_logger.ServiceLogger(log_file=log_file, error_log=Path(odir, "error.log"), log_level='INFO', log_console=True)
        print(logger)
        print(class_logger.get_instance())

    @classmethod
    def tearDownClass(cls) -> None:
        class_logger.get_instance().close()

    def test_product(self):
        """
            unitary test for product class
        """

        nc_file = TestProduct.test_path + '/../inputs/S2B_MSIL1C_20220929T103729_N0510_R008_T31TFJ_20240726T034550.nc'

        # instantiate product
        prod = Product(xr.open_dataset(nc_file))

        self.assertEqual(prod.sensor, 'S2B')
        self.assertEqual(prod.date_str, '2022-09-29T10:37:29.024Z')
        self.assertEqual(prod.width, 1229)
        self.assertEqual(prod.height, 1229)
        self.assertAlmostEqual(prod.lonmin, 4.941688508102381, 16)
        self.assertAlmostEqual(prod.lonmax, 5.253012614533291, 16)
        self.assertAlmostEqual(prod.latmin, 43.353871355515466, 16)
        self.assertAlmostEqual(prod.latmax, 43.58062066781994, 16)
        self.assertAlmostEqual(prod.xmin, 657340.0, 1)
        self.assertAlmostEqual(prod.xmax, 681920.0, 1)
        self.assertAlmostEqual(prod.ymin, 4802560.0, 1)
        self.assertAlmostEqual(prod.ymax, 4827140.0, 1)
#        self.assertAlmostEqual(float(prod.U), 1.0022856395562, 16)
        class_logger.get_instance().close()

        #load_auxiliary_data used in init, set class attributes
        # Other method not tested
        # get_flag --> used in load_flag not used
        # set_outfile --> not used
        # set_aeronetfile --> not used
        # get_elevation --> not used
        # load_flags --> not used
# class algo --> not used        
# get_elevation --> not used
