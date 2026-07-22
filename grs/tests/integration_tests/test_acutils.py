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
import numpy
from grs import acutils
from grs import class_logger
from grs.product import Product
from grs import CamsProduct
import GRSdriver


class TestAcutils(unittest.TestCase):
    """
        class for unitary test of acutils module
    """

    PROD = None
    CAMS = None
    test_path = ""

    @classmethod
    def setUpClass(cls) -> None:
        cls.test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = cls.test_path + '/../output/log_file.log'
        odir = cls.test_path + '/../output/'
        class_logger.ServiceLogger(log_file=log_file, output_dir=odir, log_level='INFO', log_console=True)

        cls.init_prod_data()
        cls.init_cams_data()

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

    @classmethod
    def init_cams_data(cls):
        cams_file = cls.test_path + '/../inputs/2023-10-12-cams-global-atmospheric-composition-forecasts.nc'
        cls.CAMS = CamsProduct(cls.PROD.raster, cams_file=cams_file)
        cls.CAMS.load()

#    def test_lut(self):
#        """
#            unitary test for lut class
#        """
#
#        # instantiate lut
#        lutf = acutils.lut(TestAcutils.PROD.band_names)
#        # execute load_lut
#        lutf.load_lut(TestAcutils.PROD.lutfine, TestAcutils.PROD.sensordata.indband)
#
#        # Verify lutf attributes
#        self.assertEqual(lutf.smac_bands, [])
#        self.assertEqual(lutf.N, 11)
#        self.assertEqual(lutf.lut_generator, 'OSOAA_h')
#        ref_wl = numpy.array([443., 490., 560., 665., 705., 740., 783., 842., 865., 1610., 2190.])
#        numpy.testing.assert_almost_equal(lutf.wl, ref_wl, 1)
#        ref_cext = numpy.array([0.08790608, 0.07532224, 0.06076635, 0.0438141, 0.03883674, 0.03482638, 0.03075605, 0.02676323, 0.02433568, 0.00419922, 0.00148824])
#        numpy.testing.assert_almost_equal(lutf.Cext, ref_cext, 8)
#        numpy.testing.assert_almost_equal(lutf.Cext550, 0.06265462594822929, 17)
#        numpy.testing.assert_array_equal(lutf.Csca, [])
#        numpy.testing.assert_array_equal(lutf.Csca550, 0)
#        ref_vza = numpy.array([0., 1.14, 2.62, 4.11, 5.61, 7.1, 8.59, 10.09, 11.58, 13.07, 14.57, 16.06, 17.55, 19.05])
#        numpy.testing.assert_almost_equal(lutf.vza, ref_vza, 2)
#        ref_sza = numpy.array([0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., 52. ,54., 56., 58., 60., 62., 64., 66., 68.])
#        numpy.testing.assert_almost_equal(lutf.sza, ref_sza, 1)
#        ref_azi = numpy.array([0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65.,
#  70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135.,
# 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 205.,
# 210., 215., 220., 225., 230., 235., 240., 245., 250., 255., 260., 265., 270., 275.,
# 280., 285., 290., 295., 300., 305., 310., 315., 320., 325., 330., 335., 340., 345.,
# 350., 355., 360.])
#        numpy.testing.assert_almost_equal(lutf.azi, ref_azi, 1)
#        ref_aot = numpy.array([0.01, 0.05, 0.1, 0.3, 0.5, 0.8])
#        numpy.testing.assert_almost_equal(lutf.aot, ref_aot, 2)
#
#        #plouf
#        # Test if log file is created
#        #self.assertTrue(os.path.isfile(log_file))


    def test_gaseous_transmittance(self):
        """
            unitary test for gaseous_transmittance class
        """
        # Instanciate gaseous_transmittance
        gaseous_transmittance_instance = acutils.GaseousTransmittance(TestAcutils.PROD, TestAcutils.CAMS)

        print(gaseous_transmittance_instance.SRF)
        ref_srf = numpy.array([0.01448181, 0.03422251, 0.07346335, 0.15444843, 0.31661424, 0.55322278,
                               0.74859405, 0.84890306, 0.89772218, 0.9215368,  0.92572844, 0.91122687,
                               0.88818926, 0.86523753, 0.84718186, 0.83875722, 0.84459078, 0.86219651,
                               0.88838714, 0.92443234, 0.96017975, 0.98685515, 1.,         0.99860078,
                               0.98076475, 0.94522089, 0.8981778,  0.85580325, 0.81841731, 0.78862047,
                               0.76460654, 0.74963742, 0.7505511,  0.76137888, 0.78244478, 0.79890084,
                               0.81016958, 0.81408888, 0.77358598, 0.62881064, 0.40397555, 0.21542098,
                               0.10715281, 0.04792877, 0.01848693, 0.00108588])
        numpy.testing.assert_almost_equal(gaseous_transmittance_instance.SRF.values[2][138:184], ref_srf, 8)

        print(gaseous_transmittance_instance.Tg_tot_coarse) 

        self.assertAlmostEqual(gaseous_transmittance_instance.xmin, 300000.0, places=1)
        self.assertAlmostEqual(gaseous_transmittance_instance.ymin, 4790220.0, places=1)
        self.assertAlmostEqual(gaseous_transmittance_instance.xmax, 409800.0, places=1)
        self.assertAlmostEqual(gaseous_transmittance_instance.ymax, 4900020.0, places=1)
        self.assertAlmostEqual(gaseous_transmittance_instance.gas_lut.wl.values[10], 351.772064, places=6)
        self.assertAlmostEqual(gaseous_transmittance_instance.gas_lut.ch4.values[20000], 0.00717464667299828, places=16)
        self.assertAlmostEqual(gaseous_transmittance_instance.Twv_lut.Twv.values[10][10][10], 0.9994355799942988, places=16)
        self.assertAlmostEqual(gaseous_transmittance_instance.air_mass_mean.values, 2.62716381, places=8)
        self.assertAlmostEqual(gaseous_transmittance_instance.pressure.values[5][5], 1000.3613857341961, places=16)
        self.assertAlmostEqual(gaseous_transmittance_instance.coef_abs_scat['h2o'], 0.3, places=1)

        # get_gaseous_transmittance
        tg_raster = gaseous_transmittance_instance.get_gaseous_transmittance()
        self.assertAlmostEqual(tg_raster.values[10][10][10], 0.002012189315168797, places=16)

        # Tgas_background
        tgas_background = gaseous_transmittance_instance.Tgas_background()
        self.assertAlmostEqual(tgas_background.values[10][10][10], 0.9999832067551963, places=16)

        # Other method not tested
        # get_gaseous_optical_thickness --> used in get_gaseous_transmittance_old so not used
        # get_gaseous_transmittance_old --> not used
        # other_gas_correction --> not used
        # water_vapor_correction --> not used
        # get_wv_transmittance_raster --> not used

    def test_misc(self):
        """
            unitary test for misc class
        """

        # test get_pressure function
        atl = 1000.0
        psl = 998.0
        palt = acutils.Misc.get_pressure(atl, psl)

        self.assertAlmostEqual(palt, 885.236756238, places=8)
