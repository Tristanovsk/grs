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
from grs import class_logger
from grs.product import Product
from grs import CamsProduct
import GRSdriver


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
        odir = cls.test_path + '/../output/'
        class_logger.ServiceLogger(log_file=log_file, output_dir=odir, log_level='INFO', log_console=True)

        cls.init_prod_data()

    @classmethod
    def tearDownClass(cls) -> None:
        class_logger.get_instance().close()

    @classmethod
    def init_prod_data(cls):
        # Load data to instanciate class
        data_file = cls.test_path + '/../inputs/S2A_MSIL1C_20231012T104951_N0509_R051_T31TCJ_20231012T143114.SAFE'
        resolution = 20
        l1c = GRSdriver.Sentinel2Driver(data_file, resolution=resolution)
        l1c.load_product()
        cls.PROD = Product(l1c.prod)

    def test_cams_products(self):
        """
            unitary test for lut class
        """
        # instantiate cams
        cams_file = '/work/datalake/watcal/ECMWF/CAMS/2023/10/12/2023-10-12-cams-global-atmospheric-composition-forecasts.nc'
        cams = CamsProduct(TestCams.PROD.raster, cams_file=cams_file)
        cams.load()

        print(cams.raster.x.values[0])
        print(cams.raster.x.values[-1])
        print(cams.raster.y.values[0])
        print(cams.raster.y.values[-1])
        self.assertAlmostEqual(cams.raster.x.values[0], 300000.0, 1) # xmin
        self.assertAlmostEqual(cams.raster.x.values[-1], 409800.0, 1) # xmax
        self.assertAlmostEqual(cams.raster.y.values[-1], 4790220.0, 1) # ymax
        self.assertAlmostEqual(cams.raster.y.values[0], 4900020.0, 1) # ymin
        
        self.assertEqual(cams.variables, ['v10', 't2m', 'msl', 'sp', 'ssa1020', 'ssa1240', 'ssa1640', 'ssa2130', 'ssa355', 'ssa380', 'ssa400', 'ssa440', 'ssa500', 'ssa550', 'ssa645', 'ssa670', 'ssa800', 'ssa865', 'aod1020', 'aod1064', 'aod1240', 'aod1640', 'aod2130', 'aod355', 'aod380', 'aod400', 'aod440', 'aod469', 'aod500', 'aod550', 'aod645', 'aod670', 'aod800', 'aod865', 'tcco', 'tchcho', 'tc_oh', 'tc_ch4', 'tcno2', 'gtco3', 'tc_c3h8', 'tcwv', 'u10'])
        print(cams.raster.u10.values[6][6])
        print(cams.raster.v10.values[6][6])
        print(cams.raster.t2m.values[6][6])
        print(cams.raster.gtco3.values[6][6])
        print(cams.raster.tcno2.values[6][6])
        #self.assertAlmostEqual(cams.cams_ssa.values[10][10][10], 0.9677208736679962, 7)
        self.assertAlmostEqual(cams.cams_aod.values[4][10][10], 0.017030397059077572, 7)
        self.assertAlmostEqual(cams.raster.u10.values[6][6], -0.993277515050355, 7)
        self.assertAlmostEqual(cams.raster.v10.values[6][6], 1.5300712706261923, 7)
        self.assertAlmostEqual(cams.raster.t2m.values[6][6], 298.96904073387566, 7)
        self.assertAlmostEqual(cams.raster.gtco3.values[6][6], 0.005542143138166818, 7)
        self.assertAlmostEqual(cams.raster.tcno2.values[6][6], 2.8778440127438184e-06, 7)


    # method not tested
    # subset_xr not used
    # download_erainterim not used
    # load_cams_data not used
    # class aeronet not used
