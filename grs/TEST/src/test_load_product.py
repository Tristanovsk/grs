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
from grs import class_logger
from grs.drivers import driver_S2_SAFE


class TestDriverS2Safe(unittest.TestCase):
    """
        class for unitary test of driver_S2_SAFE module
    """

    def init_log(self):
        test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = test_path + '/../output/log_file.log'
        odir = test_path + '/../output/'
        error_log_file = test_path + '/../output/error.log'
        # Init logger service (avoid to be spam by DEBUG log from rasterio...)
        class_logger.ServiceLogger(log_file=log_file, output_dir=odir, log_level='INFO', log_console=True)

    def test_load_product(self):
        """
            unitary test for load_product method
        """

        # Prepare log
        self.init_log()
        # Load data to instanciate class
        data_file = '/home/mp/penardc/scratch/PROJETS/OBS2CO/test/data/S2A_MSIL1C_20231012T104951_N0509_R051_T31TCJ_20231012T143114.SAFE'
        bandIds = range(13)
        resolution = 20
        # __init__
        l1c = driver_S2_SAFE.s2image(data_file, band_idx=bandIds, resolution=resolution)
        l1c.load_product()
        ref_srf = numpy.array([0.01448181, 0.03422251, 0.07346335, 0.15444843, 0.31661424, 0.55322278,
                               0.74859405, 0.84890306, 0.89772218, 0.9215368,  0.92572844, 0.91122687,
                               0.88818926, 0.86523753, 0.84718186, 0.83875722, 0.84459078, 0.86219651,
                               0.88838714, 0.92443234, 0.96017975, 0.98685515, 1.,         0.99860078,
                               0.98076475, 0.94522089, 0.8981778,  0.85580325, 0.81841731, 0.78862047,
                               0.76460654, 0.74963742, 0.7505511,  0.76137888, 0.78244478, 0.79890084,
                               0.81016958, 0.81408888, 0.77358598, 0.62881064, 0.40397555, 0.21542098,
                               0.10715281, 0.04792877, 0.01848693, 0.00108588])
        numpy.testing.assert_almost_equal(l1c.prod.SRF.values[2][138:184], ref_srf, 8)
        # Test on x and y vector
        self.assertEqual(l1c.prod.x.__len__(), 5490)
        self.assertEqual(l1c.prod.y.__len__(), 5490)
        # Test on lat/lon values
        lonmin, latmin, lonmax, latmax = l1c.prod.rio.transform_bounds(4326, recalc=True)
        xmin, ymin, xmax, ymax = l1c.prod.rio.bounds(recalc=True)
        print("lonmin:" + str(lonmin))
        print("latmin:" + str(latmin))
        print("lonmax:" + str(lonmax))
        print("latmax:" + str(latmax))
        print("xmin:" + str(xmin))
        print("ymin:" + str(ymin))
        print("xmax:" + str(xmax))
        print("ymax:" + str(ymax))
        self.assertAlmostEqual(lonmin, 0.4959285929146376, 16)
        self.assertAlmostEqual(lonmax, 1.8886830141969135, 16)
        self.assertAlmostEqual(latmin, 43.238262552870076, 16)
        self.assertAlmostEqual(latmax, 44.247830683507615, 16)
        self.assertAlmostEqual(xmin, 300000.0, 1)
        self.assertAlmostEqual(xmax, 409800.0, 1)
        self.assertAlmostEqual(ymin, 4790220.0, 1)
        self.assertAlmostEqual(ymax, 4900020.0, 1)


        # set_crs (static_method) --> not tested
