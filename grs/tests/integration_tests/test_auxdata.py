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


import unittest


class TestAuxdata(unittest.TestCase):
    """
        class for unitary test of auxdata module
    """

    def init_prod_data(self):
        # Create ref data for assert
        self.ref_S2A = dict()
        self.ref_S2A['central_wavelength'] = [442.7316, 492.4410, 559.8538, 664.6208, 704.1223,
                                              740.4838, 782.7510, 832.7700, 864.7027, 1613.6594, 2202.3662]
        self.ref_S2A['resolution'] = 20
        self.ref_S2A['lutname'] = 'S2A/lut_'
        self.ref_S2A['vza_name'] = 'view_zenith_'
        self.ref_S2A['azi_name'] = 'view_azimuth_'

        self.ref_S2B = dict()
        self.ref_S2B['central_wavelength'] = [442.3110, 492.1326, 558.9499, 664.9380, 703.8308, 739.1290,
                                              779.7236, 832.9462, 863.9796, 1610.4191, 2185.6988]
        self.ref_S2B['resolution'] = 20
        self.ref_S2B['lutname'] = 'S2B/lut_'
        self.ref_S2B['vza_name'] = 'view_zenith_'
        self.ref_S2B['azi_name'] = 'view_azimuth_'

        self.ref_LANDSAT_4 = dict()
        self.ref_LANDSAT_4['central_wavelength'] = [485.9919, 571.2153, 659.8436, 839.3312,
                                                    1679.8097, 2216.9931]
        self.ref_LANDSAT_4['resolution'] = 30
        self.ref_LANDSAT_4['lutname'] = 'L4/lut_L4_'
        self.ref_LANDSAT_4['vza_name'] = 'Zenith_'
        self.ref_LANDSAT_4['azi_name'] = 'Azimuth_'

        self.ref_LANDSAT_5 = dict()
        self.ref_LANDSAT_5['central_wavelength'] = [486.2534, 570.5804, 660.6101, 838.1482, 1677.1762, 2217.3553]
        self.ref_LANDSAT_5['resolution'] = 30
        self.ref_LANDSAT_5['lutname'] = 'L5/lut_L5_'
        self.ref_LANDSAT_5['vza_name'] = 'Zenith_'
        self.ref_LANDSAT_5['azi_name'] = 'Azimuth_'

        self.ref_LANDSAT_7 = dict()
        self.ref_LANDSAT_7['central_wavelength'] = [478.7157, 561.0342, 661.4387, 834.5691, 1650.2731, 2208.1606]
        self.ref_LANDSAT_7['resolution'] = 30
        self.ref_LANDSAT_7['lutname'] = 'L7/lut_L7_'
        self.ref_LANDSAT_7['vza_name'] = 'Zenith_'
        self.ref_LANDSAT_7['azi_name'] = 'Azimuth_'

        self.ref_LANDSAT_8 = dict()
        self.ref_LANDSAT_8['central_wavelength'] = [442.9821, 482.5889, 561.3343, 591.6667, 654.6084, 864.5711,
                                                    1609.0906, 2201.2492]
        self.ref_LANDSAT_8['resolution'] = 30
        self.ref_LANDSAT_8['lutname'] = 'L8/lut_L8_'
        self.ref_LANDSAT_8['vza_name'] = 'Zenith_'
        self.ref_LANDSAT_8['azi_name'] = 'Azimuth_'

    def verif_values(self, sensordata, ref):
        # Only a sample of values is tested
        self.assertAlmostEqual(sensordata.central_wavelength, ref['central_wavelength'], 4)
        self.assertEqual(sensordata.resolution, ref['resolution'])
        self.assertEqual(sensordata.lutname, ref['lutname'])
        self.assertEqual(sensordata.vza_name, ref['vza_name'])
        self.assertEqual(sensordata.azi_name, ref['azi_name'])

    # def test_sensordata(self):
    #     """
    #         Obsolete, need to be removed ( SensorData from auxdata is no longer a valid code)
    #         unitary test for sensordata class
    #     """
    #
    #     # Prepare test data
    #     self.init_prod_data()
    #     # instantiate sensordata
    #     sensor='S2A'
    #     sensordata = auxdata.SensorData(sensor)
    #     self.verif_values(sensordata, self.ref_S2A)
    #     sensor='S2B'
    #     sensordata = auxdata.SensorData(sensor)
    #     self.verif_values(sensordata, self.ref_S2B)
    #     sensor='LANDSAT_4'
    #     sensordata = auxdata.SensorData(sensor)
    #     self.verif_values(sensordata, self.ref_LANDSAT_4)
    #     sensor='LANDSAT_5'
    #     sensordata = auxdata.SensorData(sensor)
    #     self.verif_values(sensordata, self.ref_LANDSAT_5)
    #     sensor='LANDSAT_7'
    #     sensordata = auxdata.SensorData(sensor)
    #     self.verif_values(sensordata, self.ref_LANDSAT_7)
    #     sensor='LANDSAT_8'
    #     sensordata = auxdata.SensorData(sensor)
    #     self.verif_values(sensordata, self.ref_LANDSAT_8)

# Class not used
#    def test_Aeronet(self):
#        """
#            unitary test for Aeronet class
#        """
