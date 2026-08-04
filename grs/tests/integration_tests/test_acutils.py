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
from pathlib import Path
from grs import acutils
from grs import class_logger
from grs.product import Product
from grs import CamsProduct
import os
import xarray as xr

class TestAcutils(unittest.TestCase):
    """
        class for unitary test of acutils module
    """

    PROD = None
    CAMS = None
    test_path = ""


    @classmethod
    def setUpClass(cls) -> None:
        print("setUpClass", flush=True)
        cls.test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = cls.test_path + '/../output/log_file.log'
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        odir = cls.test_path + '/../output/'
        class_logger.ServiceLogger(log_file=log_file, error_log=Path(odir, "error.log"), log_level='INFO', log_console=True)

        cls.init_prod_data_nc()
        cls.init_cams_data()
        print("END setUpClass", flush=True)

    @classmethod
    def tearDownClass(cls) -> None:
        class_logger.get_instance().close()

    @classmethod
    def init_prod_data_nc(cls):
        nc_file =  cls.test_path + '/../inputs/S2B_MSIL1C_20220929T103729_N0510_R008_T31TFJ_20240726T034550.nc'
        cls.PROD = Product(xr.open_dataset(nc_file))
        

    @classmethod
    def init_cams_data(cls):
        cams_file = cls.test_path + '/../inputs/cams_forecast_2022-09-29.nc'
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
        print("test_gaseous_transmittance", flush=True)
        # Instanciate gaseous_transmittance
        gaseous_transmittance_instance = acutils.GaseousTransmittance(TestAcutils.PROD, TestAcutils.CAMS)
        print("After construction de acutils.GaseousTransmittance", flush=True)
        print(gaseous_transmittance_instance.SRF)
        ref_srf = numpy.array([0.03660414, 0.08100583, 0.16917887, 0.33278275, 0.58622795, 0.8091641,
                                0.913051, 0.94472283, 0.94898814, 0.9436913, 0.9284567, 0.9125694,
                                0.9007804, 0.89958596, 0.9054714, 0.92045355, 0.94065666, 0.9619968,
                                0.98186743, 0.9985841, 1., 0.99279886, 0.9780133, 0.95301175,
                                0.9266333, 0.8935913, 0.8694179, 0.84827, 0.839083, 0.83206207,
                                0.8291787, 0.8330584, 0.84630936, 0.86396307, 0.8726808, 0.8681834,
                                0.8554947, 0.80839056, 0.6765088, 0.45584205, 0.24737576, 0.12765466,
                                0.0589016, 0.02564742, 0.00515905, numpy.nan])

        print(ref_srf)
        print(gaseous_transmittance_instance.SRF.values[2][138:184])

        diff = numpy.abs(gaseous_transmittance_instance.SRF.values[2][138:184] - ref_srf)
        print("Différence max :", numpy.nanmax(diff))

        idx = numpy.where(diff > 0)[0]
        print("Indices différents :", idx)

        for i in idx:
            print(i, repr(gaseous_transmittance_instance.SRF.values[2][138:184][i]), repr(ref_srf[i]), 
                  repr(gaseous_transmittance_instance.SRF.values[2][138:184][i] - ref_srf[i]))

        numpy.testing.assert_almost_equal(gaseous_transmittance_instance.SRF.values[2][138:184], ref_srf, 8)

        print(gaseous_transmittance_instance.Tg_tot_coarse) 

        print("xmin=",gaseous_transmittance_instance.xmin)
        print("ymin=",gaseous_transmittance_instance.ymin)
        print("xmax=",gaseous_transmittance_instance.xmax)
        print("ymax=",gaseous_transmittance_instance.ymax)
        print("gas_lut.wl.values[10]=",gaseous_transmittance_instance.gas_lut.wl.values[10])
        print("gas_lut.wl.ch4.values[20000]=",gaseous_transmittance_instance.gas_lut.ch4.values[20000])
        print("gas_lut.wl.Twv.values[10][10][10]=",gaseous_transmittance_instance.gas_lut.Twv.values[10][10][10])
        print("air_mass_mean.values=",gaseous_transmittance_instance.air_mass_mean.values)
        print("pressure.values[5][5]=",gaseous_transmittance_instance.pressure.values[5][5])
        print("pressure.coef_abs_scat['h2o']=",gaseous_transmittance_instance.pressure.coef_abs_scat['h2o'])

        print("")
        print("tg_raster.values[10][10][10]=", tg_raster.values[10][10][10])
        print("tgas_background.values[10][10][10]=", tgas_background.values[10][10][10])
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
        print("test_misc", flush=True)
        palt = acutils.Misc.get_pressure(atl, psl)
        print("After acutils.Misc.get_pressure", flush=True)
        self.assertAlmostEqual(palt, 885.236756238, places=8)

