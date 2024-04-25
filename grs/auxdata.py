'''
Module dedicated to auxiliary data management.
'''

import os
import logging
import numpy as np
import pandas as pd
import xarray as xr

from scipy.interpolate import interp1d

from pkg_resources import resource_filename
from importlib.resources import files

opj = os.path.join

# ------------------------
# set threshold for masking
O2band_cloud = [0.08, 0.12]
high_nir_threshold = 0.0275
hcld_threshold = 0.003  # (0.2 % check http://www.cesbio.ups-tlse.fr/multitemp/?p=12894)

# -----------------
# set values bracketing the Normalized Difference Water Index for rough land/water masking
NDWI_nir_threshold = [-0.03, 1.]
NDWI_swir_threshold = [0.12, 2.]

# NDWI_nir_threshold = [-0., 1.]


# ******************************************************************************************************
dir, filename = os.path.split(__file__)

sunglint_eps_file = files('grs.data.aux').joinpath('mean_rglint_small_angles_vza_le_12_sza_le_60.txt')
rayleigh_file = files('grs.data.aux').joinpath('rayleigh_bodhaine.txt')


class AuxData():
    '''
    Class to load auxiliary data:

    - Rayleigh optical thickness
    - Factor for the spectral variation of the Fresnel reflectance to be applied to sunglint calculation.
    '''

    def __init__(self, wl=None):
        # load data from raw files

        self.sunglint_eps = pd.read_csv(sunglint_eps_file, sep='\s+', index_col=0).to_xarray()
        self.rayleigh()

        # reproject onto desired wavelengths
        if wl is not None:
            self.sunglint_eps = self.sunglint_eps['mean'].interp(wl=wl)
            self.rot = self.rot.interp(wl=wl)

    def rayleigh(self):
        '''
        Rayleigh Optical Thickness for
        P=1013.25mb,
        T=288.15K,
        CO2=360ppm
        from
        Bodhaine, B.A., Wood, N.B, Dutton, E.G., Slusser, J.R. (1999). On Rayleigh
        Optical Depth Calculations, J. Atmos. Ocean Tech., 16, 1854-1861.
        '''
        data = pd.read_csv(rayleigh_file, skiprows=16, sep=' ', header=None)
        data.columns = ('wl', 'rot', 'dpol')
        self.rot_pressure_ref = 101325
        self.rot = data.set_index('wl').to_xarray().rot


class SensorData:
    ''' TODO Obsolete module: needs to be removed

    Dictionnaries of the auxilliary data, functions to be applied for the sensors data to process.

    Arguments:
        * ``sensor`` -- Name of the sensor to process:
            * LANDSAT_4
            * LANDSAT_5
            * LANDSAT_7
            * LANDSAT_8
            * S2A
            * S2B
        * ``indband`` -- indexes of the requested bands (bands that will be processed)

    Construct object with values:
        * ``rot`` -- Rayleigh optical thickness for each band
        * ``band names`` for angles files format
        * ``rg`` -- Cox-Munk Fresnel reflection factor ratio (R(wl)/Rref(2190nm) [Harmel et al., 2018]
        * ``angle_generator`` -- function to compute angles from metadata
    '''

    def __init__(self, sensor):
        from . import config as cfg

        self.sensor = sensor

        ##################################
        # Sensors characteristics
        ##################################
        # Computed from Thuillier 2003 data set [Thuillier, G., M. Herse, P. C. Simon, D. Labs, H. Mandel, D. Gillotay,
        #               and T. Foujols, 2003, 'The solar spectral irradiance from 200
        #               to 2400 nm as measured by the SOLSPEC spectrometer from the
        #               ATLAS 1-2-3 and EURECA missions, Solar Physics, 214(1): 1-22
        #    [u'CoastalAerosol', u'Blue', u'Green', u'Red', u'NIR', u'Cirrus', u'SWIR1', u'SWIR2', u'Pan']
        #    [1895.42490193, 2004.57767987, 1820.64535207, 1549.46770711, 951.72893589, 366.97254025,
        #     247.55307624, 85.46266307, 1723.80908067]
        LANDSAT8_ESUN = [1895.42490193, 2004.57767987, 1820.64535207, 1723.80908067, 1549.46770711,
                         951.72893589, 247.55307624, 85.46266307]

        # From Chander, G., Markham, B.L., Helder, D.L. (2008)
        # Summary of current radiometric calibration coefficients for Landsat MSS, TM, ETM+, and EO-1 ALI sensors
        # http://landsathandbook.gsfc.nasa.gov/pdfs/Landsat_Calibration_Summary_RSE.pdf
        LANDSAT4_ESUN = [1983.0, 1795.0, 1539.0, 1028.0, 219.8, 83.49]
        LANDSAT5_ESUN = [1983.0, 1796.0, 1536.0, 1031.0, 220.0, 83.44]
        LANDSAT7_ESUN = [1997.0, 1812.0, 1533.0, 1039.0, 230.8, 84.90]

        InfoSat = {
            'S2A': {'sensor': 'S2A',
                    'resolution': 20,
                    'indband': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    'indlut': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    'central_wavelength': [442.7316, 492.4410, 559.8538, 664.6208, 704.1223,
                                           740.4838, 782.7510, 832.7700, 864.7027, 1613.6594, 2202.3662],
                    'solar_irr': None,
                    'band_names': ('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'),
                    'ang_names': ('B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12'),
                    'vza_name': 'view_zenith_',
                    'azi_name': 'view_azimuth_',
                    'ang_proc': None,
                    'cirrus': ['B10', hcld_threshold],
                    'lut_name': 'S2A/lut_',
                    # rot and rglint from old RSR (now updated since Jan-2018)
                    # 'rot': [0.23385, 0.14944, 0.09034, 0.0449976, 0.035565,
                    #         0.02901, 0.023186, 0.018123,0.0154835, 0.0012656, 0.000365],
                    # 'rglint': [1.286189, 1.266816, 1.249609, 1.230385, 1.224821,
                    #            1.220275, 1.215545, 1.209850, 1.206606, 1.124581, 1.000000],
                    'rot': [0.23650534, 0.15478519, 0.09043935, 0.04495076, 0.0355168,
                            0.02896814, 0.02315269, 0.01830901, 0.01549064, 0.00126556, 0.00036496],
                    'rglint': [1.28672512, 1.26815692, 1.24964349, 1.23035937, 1.22478961,
                               1.2202435, 1.21551569, 1.21009995, 1.206616, 1.12458056, 1.],
                    'ndwi_nir_conf': [2, 8, NDWI_nir_threshold],
                    'ndwi_swir_conf': [8, 9, NDWI_swir_threshold],
                    'high_nir': [8, high_nir_threshold],
                    'O2band': ['B9', *O2band_cloud],
                    'l1_flags': ['opaque_clouds_10m', 'cirrus_clouds_10m', '']
                    },

            # TODO update for S2B: lut_name, rot, rglint for RSR
            'S2B': {'sensor': 'S2B',
                    'resolution': 20,
                    'indband': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    'indlut': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    'central_wavelength': [442.3110, 492.1326, 558.9499, 664.9380, 703.8308, 739.1290,
                                           779.7236, 832.9462, 863.9796, 1610.4191, 2185.6988],
                    'solar_irr': None,
                    'band_names': ('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'),
                    'ang_names': ('B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12'),
                    'vza_name': 'view_zenith_',
                    'azi_name': 'view_azimuth_',
                    'ang_proc': None,
                    'cirrus': ['B10', hcld_threshold],
                    'lut_name': 'S2B/lut_',
                    'rot': [0.23745233, 0.15521662, 0.09104522, 0.04485828, 0.03557663,
                            0.02918357, 0.02351827, 0.01829234, 0.01554304, 0.00127606, 0.00037652],
                    'rglint': [1.27878376, 1.26024821, 1.24195311, 1.22253644, 1.21708874, 1.21269296,
                               1.20815455, 1.20243389, 1.19906708, 1.11793355, 1.],
                    'ndwi_nir_conf': [2, 8, NDWI_nir_threshold],
                    'ndwi_swir_conf': [8, 9, NDWI_swir_threshold],
                    'high_nir': [8, high_nir_threshold],
                    'O2band': ['B9', *O2band_cloud],
                    'l1_flags': ['opaque_clouds_10m', 'cirrus_clouds_10m', '']
                    },

            'LANDSAT_4': {'sensor': 'LANDSAT_4',
                          'resolution': 30,
                          'indband': [0, 1, 2, 3, 4, 5],
                          'indlut': [0, 1, 2, 3, 4, 5],
                          'central_wavelength': [485.9919, 571.2153, 659.8436, 839.3312,
                                                 1679.8097, 2216.9931],
                          'solar_irr': LANDSAT4_ESUN,
                          'band_names': (
                              'radiance_1', 'radiance_2', 'radiance_3', 'radiance_4', 'radiance_5', 'radiance_7'),
                          'ang_names': ('B01', 'B02', 'B03', 'B04', 'B05', 'B07'),
                          'vza_name': 'Zenith_',
                          'azi_name': 'Azimuth_',
                          'ang_proc': os.path.join(cfg.grs_root, 'landsat_angles/TM/landsat_angles'),
                          'cirrus': ['no', hcld_threshold],
                          'lut_name': 'L4/lut_L4_',
                          'rot': [0.16389, 0.084813, 0.046737, 0.01785, 0.001091, 0.0003578],
                          'rglint': [1.28037, 1.25722, 1.24089, 1.21893, 1.12406, 1.0],
                          'ndwi_nir_conf': [1, 3, NDWI_nir_threshold],
                          'ndwi_swir_conf': [3, 4, NDWI_swir_threshold],
                          'high_nir': [3, high_nir_threshold],
                          'O2band': None,
                          'l1_flags': ['', '', '']
                          },

            'LANDSAT_5': {'sensor': 'LANDSAT_5',
                          'resolution': 30,
                          'indband': [0, 1, 2, 3, 4, 5],
                          'indlut': [0, 1, 2, 3, 4, 5],
                          'central_wavelength': [486.2534, 570.5804, 660.6101, 838.1482, 1677.1762, 2217.3553],
                          'solar_irr': LANDSAT5_ESUN,
                          'band_names': (
                              'radiance_1', 'radiance_2', 'radiance_3', 'radiance_4', 'radiance_5', 'radiance_7'),
                          'ang_names': ('B01', 'B02', 'B03', 'B04', 'B05', 'B07'),
                          'vza_name': 'Zenith_',
                          'azi_name': 'Azimuth_',
                          'ang_proc': os.path.join(cfg.grs_root, 'landsat_angles/TM/landsat_angles'),
                          'cirrus': ['no', hcld_threshold],
                          'lut_name': 'L5/lut_L5_',
                          'rot': [0.163419, 0.085201, 0.046514, 0.017945, 0.001098, 0.0003576],
                          'rglint': [1.28046, 1.257556, 1.24096, 1.21924, 1.12462, 1.0],
                          'ndwi_nir_conf': [1, 3, NDWI_nir_threshold],
                          'ndwi_swir_conf': [3, 5, NDWI_swir_threshold],
                          'high_nir': [3, high_nir_threshold],
                          'O2band': None,
                          'l1_flags': ['', '', '']
                          },

            'LANDSAT_7': {'sensor': 'LANDSAT_7',
                          'resolution': 30,
                          'indband': [0, 1, 2, 3, 4, 5],
                          'indlut': [0, 1, 2, 3, 4, 5],
                          'central_wavelength': [478.7157, 561.0342, 661.4387, 834.5691, 1650.2731, 2208.1606],
                          'solar_irr': LANDSAT7_ESUN,
                          'band_names': (
                              'radiance_1', 'radiance_2', 'radiance_3', 'radiance_4', 'radiance_5', 'radiance_7'),
                          'ang_names': ('B01', 'B02', 'B03', 'B04', 'B05', 'B07'),
                          'vza_name': 'Zenith_',
                          'azi_name': 'Azimuth_',
                          'ang_proc': os.path.join(cfg.grs_root, 'landsat_angles/TM/landsat_angles'),
                          'cirrus': ['no', hcld_threshold],
                          'lut_name': 'L7/lut_L7_',
                          'rot': [0.174702, 0.091113, 0.046100, 0.018231, 0.001169, 0.0003645],
                          'rglint': [1.27891, 1.25551, 1.236695, 1.215597, 1.12477, 1.0],
                          'ndwi_nir_conf': [1, 3, NDWI_nir_threshold],
                          'ndwi_swir_conf': [3, 5, NDWI_swir_threshold],
                          'high_nir': [3, high_nir_threshold],
                          'O2band': None,
                          'l1_flags': ['', '', '']
                          },

            'LANDSAT_8': {'sensor': 'LANDSAT_8',
                          'resolution': 30,
                          'indband': [0, 1, 2, 3, 4, 5, 6, 7],
                          'indlut': [0, 1, 2, 3, 4, 5, 6, 7],
                          'central_wavelength': [442.9821, 482.5889, 561.3343, 591.6667, 654.6084,
                                                 864.5711, 1609.0906, 2201.2492],
                          'solar_irr': LANDSAT8_ESUN,
                          'band_names': ('coastal_aerosol', 'blue', 'green', 'panchromatic', 'red', 'near_infrared',
                                         'swir_1', 'swir_2'),
                          'ang_names': ('B01', 'B02', 'B03', 'B09', 'B04', 'B05', 'B06', 'B07', 'B08'),
                          'vza_name': 'Zenith_',
                          'azi_name': 'Azimuth_',
                          'ang_proc': os.path.join(cfg.grs_root, 'landsat_angles/OLI/l8_angles'),
                          'cirrus': ['cirrus', hcld_threshold],
                          'lut_name': 'L8/lut_L8_',
                          'rot': [0.22899, 0.16786, 0.08999, 0.077445, 0.04785, 0.015507, 0.00128, 0.000366],
                          'rglint': [1.28637, 1.271187, 1.249180, 1.243706, 1.231680, 1.206419, 1.125028, 1.000000],
                          'ndwi_nir_conf': [2, 5, NDWI_nir_threshold],
                          'ndwi_swir_conf': [5, 6, NDWI_swir_threshold],
                          'high_nir': [5, high_nir_threshold],
                          'O2band': None,
                          'l1_flags': ['cloud_confidence_high', 'cirrus_confidence_high',
                                       'cloud_shadow_confidence_high']
                          }
        }

        info = InfoSat[sensor]
        idx = info['indband']
        # set and reorganize data for the requested bands
        self.indband = idx
        self.central_wavelength = info['central_wavelength']
        self.resolution = info['resolution']
        self.rot = np.array(info['rot'])[idx]
        self.rg = np.array(info['rglint'])[idx]
        self.lutname = info['lut_name']
        self.angle_processor = info['ang_proc']
        self.angle_names = np.array(info['ang_names'])[idx]
        self.vza_name = info['vza_name']
        self.azi_name = info['azi_name']
        self.band_names = np.array(info['band_names'])[idx]
        self.solar_irr = info['solar_irr']
        self.NDWI_vis, self.NDWI_nir, self.NDWI_threshold = info['ndwi_nir_conf']
        self.NDWI_swir_nir, self.NDWI_swir_swir, self.NDWI_swir_threshold = info['ndwi_swir_conf']
        self.cloud_flag, self.cirrus_flag, self.shadow_flag = info['l1_flags']
        self.cirrus = info['cirrus']
        self.high_nir = info['high_nir']
        self.O2band = info['O2band']
