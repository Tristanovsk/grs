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


import os
import os.path
import unittest
from pathlib import Path
from grs import class_logger
from grs.product import Product
from grs import CamsProduct
import xarray as xr


class TestCams(unittest.TestCase):
    """
        class for unitary test of cams module
    """

    PROD = None
    test_path = ""

    @classmethod
    def setUpClass(cls) -> None:
        cls.test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = cls.test_path + '/../output/log_file.log'
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        odir = cls.test_path + '/../output/'
        class_logger.ServiceLogger(log_file=log_file, error_log=Path(odir, "error.log"), log_level='INFO', log_console=True)

        cls.init_prod_data()

    @classmethod
    def tearDownClass(cls) -> None:
        class_logger.get_instance().close()

    @classmethod
    def init_prod_data(cls):
        # Load data to instanciate class
        nc_file =  cls.test_path + '/../inputs/S2B_MSIL1C_20220929T103729_N0510_R008_T31TFJ_20240726T034550.nc'
        cls.PROD = Product(xr.open_dataset(nc_file))

    def test_cams_products(self):
        """
            unitary test for lut class
        """
        # instantiate cams
        cams_file = TestCams.test_path + '/../inputs/cams_forecast_2022-09-29.nc'
        cams = CamsProduct(TestCams.PROD.raster, cams_file=cams_file)
        cams.load()

        print(cams.raster.x.values[0])
        print(cams.raster.x.values[-1])
        print(cams.raster.y.values[0])
        print(cams.raster.y.values[-1])
        self.assertAlmostEqual(cams.raster.x.values[0], 657340.0, 1) # xmin
        self.assertAlmostEqual(cams.raster.x.values[-1], 681920.0, 1) # xmax
        self.assertAlmostEqual(cams.raster.y.values[0], 4827140.0, 1) # ymin
        self.assertAlmostEqual(cams.raster.y.values[-1], 4802560.0, 1) # ymax
        
        self.assertEqual(cams.variables, ['v10', 't2m', 'msl', 'sp', 'amaod550', 'bcaod550', 'duaod550', 'niaod550', 'omaod550', 'ssaod550', 'suaod550', 'aod1240', 'aod469', 'aod550', 'aod670', 'aod865', 'tcco', 'tc_ch4', 'tcno2', 'gtco3', 'tcwv', 'u10'])
        print(cams.raster.u10.values[6][6])
        print(cams.raster.v10.values[6][6])
        print(cams.raster.t2m.values[6][6])
        print(cams.raster.gtco3.values[6][6])
        print(cams.raster.tcno2.values[6][6])

        self.assertAlmostEqual(cams.cams_aod.values[4][10][10], 0.013252163430922682, 7)
        self.assertAlmostEqual(cams.raster.u10.values[6][6], 5.255191621809645, 7)
        self.assertAlmostEqual(cams.raster.v10.values[6][6], -1.4699149996404888, 7)
        self.assertAlmostEqual(cams.raster.t2m.values[6][6], 290.8590535077976, 7)
        self.assertAlmostEqual(cams.raster.gtco3.values[6][6], 0.006902521270544851, 7)
        self.assertAlmostEqual(cams.raster.tcno2.values[6][6], 2.464762377680245e-06, 7)


    # method not tested
    # subset_xr not used
    # download_erainterim not used
    # load_cams_data not used
    # class aeronet not used
