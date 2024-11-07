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
import xarray
import datetime
import xml.etree.ElementTree as ET
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

    def test_driver_S2_SAFE(self):
        """
            unitary test for driver_S2_SAFE class
        """

        # Prepare log
        self.init_log()
        # Load data to instanciate class
        data_file = '/home/mp/penardc/scratch/PROJETS/OBS2CO/test/data/S2A_MSIL1C_20231012T104951_N0509_R051_T31TCJ_20231012T143114.SAFE'
        bandIds = range(13)
        resolution = 20
        # __init__
        l1c = driver_S2_SAFE.s2image(data_file, band_idx=bandIds, resolution=resolution)

        self.assertEqual(l1c.epsg, 32631)
        self.assertEqual(l1c.datetime, datetime.datetime(2023, 10, 12, 10, 49, 51))

        # parse_angular_grid_node (static_method)
        test_path = os.path.dirname(os.path.abspath(__file__))
        xml_granule = test_path + '/../data/MTD_TL.xml'
        with open(xml_granule) as xml_file:
            tree = ET.parse(xml_file)
            root = tree.getroot()
        raw_sza = l1c.parse_angular_grid_node(root.find('.//Tile_Angles/Sun_Angles_Grid/Zenith'))
        ref_6_raw_sza = numpy.array([30.4342, 30.4184, 30.4027, 30.3871, 30.3714, 30.3559, 30.3404, 30.3249, 30.3095, 30.2941,
         30.2788, 30.2636, 30.2484, 30.2332, 30.2181, 30.2031, 30.1881, 30.1731, 30.1582, 30.1434, 30.1286, 30.1138, 30.0992])
        numpy.testing.assert_almost_equal(raw_sza[6], ref_6_raw_sza, 4)

        # linfit (static_method)
        beta = [1, 2, 3]
        x = [0, 1]
        retour = l1c.linfit(beta, x)
        self.assertEqual(retour, 5)

        # lin2D (static_method)
        betas = numpy.array([[-3.30291994e-04,  9.61178671e-05, -2.87376467e+02], [ 7.57478910e-04, -2.16543919e-04,  9.55397893e+02], [-2.01299717e-03,  5.80951571e-04, -2.10819830e+03], [ 1.94844747e-03, -5.57771187e-04,  2.23379337e+03], [-6.61004247e-04,  1.91296643e-04, -3.62036549e+02], [ 3.37842682e-04, -9.43647195e-05,  5.85801498e+02], [-1.52735093e-04,  4.59907326e-05,  1.41830196e+02], [-1.52735093e-04,  4.59907326e-05,  1.41830196e+02]])
        x = [300000., 300020.00364365, 300040.0072873, 409759.9927127, 409779.99635635, 409800.]
        y = [4900020., 4899999.99635635, 4899979.9927127, 4790260.0072873, 4790240.00364365, 4790220.]
        mask = numpy.array([[ 4,  4,  4, 9,  9,  9], [ 4,  4,  4, 9,  9,  9], [ 4,  4,  4, 9,  9,  9], [ 5,  5,  5, 10, 10, 10], [ 5,  5,  5, 10, 10, 10], [ 5,  5,  5, 10, 10, 10]])
        new_arr = numpy.full((6, 6), numpy.nan, dtype=numpy.float32)
        l1c.lin2D(new_arr, x, y, mask, betas, detector_offset=4, scale_factor=1)
        detector_offset = 4

        beta = betas[mask - detector_offset]
        x = numpy.array(x)
        y = numpy.array(y)
        ref_arr = ((beta[:,:,0]*x)+(beta[:,:,1].transpose()*y).transpose()+beta[:,:,2])
        ref_arr = ref_arr.astype(numpy.float32, copy=False)
        numpy.testing.assert_almost_equal(new_arr, ref_arr, 5)

        # scat_angle (static_method)
        sza = 120
        vza = 80
        azi = 60
        angle = l1c.scat_angle(sza, vza, azi)
        numpy.testing.assert_almost_equal(angle, 109.85312574203562, 16)

        # data_fitting
        x0 = [300000., 304990.90909091]
        y0 = [4900020., 4895029.09090909, 4890038.18181818, 4885047.27272727, 4880056.36363636, 4875065.45454545]
        arr = numpy.array([[4.79167, 4.42863], [4.68747, 4.32531], [4.58342, numpy.nan], [4.4796, numpy.nan],
                           [4.37569, numpy.nan], [4.27251, numpy.nan]])
        data_arr = xarray.DataArray(arr)

        resfit = l1c.data_fitting(x0, y0, data_arr)
        ref = numpy.array([-7.26198111e-05, 2.08060512e-05, -7.53727940e+01])
        numpy.testing.assert_almost_equal(resfit, ref, 6)

        # get_detector_mask
        mask = l1c.get_detector_mask(bandId=0, resolution=20, detector_mask_name='DETFOO')
        self.assertEqual(mask[0,0],4) 
        self.assertEqual(mask[0,400],5) 
        self.assertEqual(mask[0,1800],6) 
        self.assertEqual(mask[0,5400],9) 

        # get_raw_angles
        ref_rsa_sza = numpy.array([52.1954, 52.1858, 52.1762, 52.1666, 52.157, 52.1474, 52.1379, 52.1284, 52.1189,
                               52.1095, 52.1001, 52.0906, 52.0813, 52.0719, 52.0626, 52.0533, 52.044,  52.0347,
                               52.0255, 52.0163, 52.0071, 51.9979, 51.9888])
        ref_rsa_sazi = numpy.array([165.662, 165.738, 165.814, 165.89, 165.966, 166.042, 166.118, 166.194, 166.27,
                                166.347, 166.423, 166.499, 166.575, 166.651, 166.728, 166.804, 166.88, 166.957,
                                167.033, 167.109, 167.186, 167.262, 167.338])
        ref_rva_vza = numpy.array([2.12028, 1.7997, 1.51169, 1.27855])
        ref_rva_vazi = numpy.array([72.7826, 66.3402, 57.2675, 44.4255])
        l1c.get_raw_angles()
        numpy.testing.assert_almost_equal(l1c.raw_sun_ang.sza.values[10][:], ref_rsa_sza, 4)
        numpy.testing.assert_almost_equal(l1c.raw_sun_ang.sazi.values[10][:], ref_rsa_sazi, 4)
        numpy.testing.assert_almost_equal(l1c.raw_view_ang.vza.values[4][2][10][4:8], ref_rva_vza, 5)
        numpy.testing.assert_almost_equal(l1c.raw_view_ang.vazi.values[4][2][10][4:8], ref_rva_vazi, 5)
        self.assertEqual(l1c.detector_num, 7) 

        # get_band_angle_as_numpy
        raw_vza = l1c.raw_view_ang.vza
        new_arr = l1c.get_band_angle_as_numpy(raw_vza, bandId=0, resolution=20)
        ref_arr = numpy.array([4.7913294, 4.789877, 4.788424])
        numpy.testing.assert_almost_equal(new_arr[0][0:3], ref_arr, 7)

        # load_bands
        l1c.load_bands(add_time=False)
        prod = l1c.prod.bands.values[2][10][100:110]
        ref_prod = numpy.array([0.100875, 0.1045, 0.105525, 0.0979, 0.11225, 0.1621, 0.17015, 0.1693, 0.15565, 0.1481])
        numpy.testing.assert_almost_equal(prod, ref_prod, 7)

        # get_all_band_angles
        l1c.get_all_band_angles(method='linear')
        ref_vza = numpy.array([4.5254145, 4.5238776, 4.5223403, 4.5208035])
        ref_raa = numpy.array([289.07602951, 289.07271129, 289.06938544, 289.06606722])
        ref_sza = numpy.array([52.63796653, 52.63792806, 52.63788958, 52.6378511 ])
        numpy.testing.assert_almost_equal(l1c.prod['vza'].values[2][10][4:8], ref_vza, 7)
        numpy.testing.assert_almost_equal(l1c.prod['raa'].values[2][10][4:8], ref_raa, 8)
        numpy.testing.assert_almost_equal(l1c.prod['sza'].values[10][4:8], ref_sza, 8)

        # load_geom call get_raw_angles and get_all_band_angles (only)
        l1c.load_geom()
        # Same assert as above, to check there is no pb with load_geom
        numpy.testing.assert_almost_equal(l1c.prod['vza'].values[2][10][4:8], ref_vza, 7)
        numpy.testing.assert_almost_equal(l1c.prod['raa'].values[2][10][4:8], ref_raa, 8)
        numpy.testing.assert_almost_equal(l1c.prod['sza'].values[10][4:8], ref_sza, 8)
        numpy.testing.assert_almost_equal(l1c.raw_sun_ang.sza.values[10][:], ref_rsa_sza, 4)
        numpy.testing.assert_almost_equal(l1c.raw_sun_ang.sazi.values[10][:], ref_rsa_sazi, 4)
        numpy.testing.assert_almost_equal(l1c.raw_view_ang.vza.values[4][2][10][4:8], ref_rva_vza, 5)
        numpy.testing.assert_almost_equal(l1c.raw_view_ang.vazi.values[4][2][10][4:8], ref_rva_vazi, 5)
        self.assertEqual(l1c.detector_num, 7) 

        # load_product : call load_bands and load_geom (+ other things)
        # delete l1c to avoid pb with load_bands method
        del l1c

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
        self.assertAlmostEqual(lonmin, 0.4959285929146376, 16)
        self.assertAlmostEqual(lonmax, 1.8886830141969135, 16)
        self.assertAlmostEqual(latmin, 43.238262552870076, 16)
        self.assertAlmostEqual(latmax, 44.247830683507615, 16)
        xmin, ymin, xmax, ymax = l1c.prod.rio.bounds(recalc=True)
        self.assertAlmostEqual(xmin, 300000.0, 1)
        self.assertAlmostEqual(xmax, 409800.0, 1)
        self.assertAlmostEqual(ymin, 4790220.0, 1)
        self.assertAlmostEqual(ymax, 4900020.0, 1)

        del l1c

        # set_crs (static_method) --> not tested
