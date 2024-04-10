'''
Module to build GRS object from input image.
'''

import os

import numpy as np
import xarray as xr
import datetime

import logging
from importlib.resources import files
import yaml
from . import AuxData, __version__, __package__

opj = os.path.join

configfile = files(__package__) / 'config.yml'
with open(configfile, 'r') as file:
    config = yaml.safe_load(file)


class Product():
    '''
    Assign the GRS product object and attributes from L1C image raster

    '''

    def __init__(self, l1c_obj=None,
                 auxdatabase='cams-global-atmospheric-composition-forecasts',
                 output='Rrs'):
        '''
        Set the metadata of GRS object.

        :param l1c_obj: input image
        :param auxdatabase: Deprecated
        :param output: Deprecated
        '''

        self.processor = __package__ + '_' + __version__

        self.raster = l1c_obj  # .prod
        # TODO check why drivers sends an object for wl coordinates instead of array of int
        self.raster['wl'] = self.raster['wl'].astype(int)
        self.sensor = self.raster.attrs['satellite']
        self.date_str = self.raster.attrs['acquisition_date']
        if 'S2' in self.sensor:
            self.date = datetime.datetime.strptime(self.date_str, '%Y-%m-%dT%H:%M:%S.%fZ')
        else:
            self.date = datetime.datetime.strptime(self.date_str, '%Y-%m-%d %H:%M:%S')

        if not 'time' in self.raster.coords:
            self.raster = self.raster.assign_coords({'time': self.date})

        # add metadata for future export to L2product
        self.raster.attrs['processing_time'] = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f')
        self.raster.attrs['processor'] = self.processor
        self.raster.attrs['version'] = __version__

        self.x = self.raster.x
        self.y = self.raster.y
        self.width = self.x.__len__()
        self.height = self.y.__len__()

        self.lonmin, self.latmin, self.lonmax, self.latmax = self.raster.rio.transform_bounds(4326)
        # set longitude between 0 and 360 deg
        self.lonmin, self.lonmax, = self.lonmin % 360, self.lonmax % 360
        self.xmin, self.ymin, self.xmax, self.ymax = self.raster.rio.bounds()

        self.wl = self.raster.wl
        self.central_wavelength = self.raster.wl_true.values

        # correct for bug with VZA == Inf
        self.raster['vza'] = self.raster.vza.where(self.raster.vza < 88)
        self.vza_mean = self.raster.vza.mean()
        self.sza_mean = self.raster.sza.mean()
        self.air_mass_mean = 1. / np.cos(np.radians(self.sza_mean)) + 1. / np.cos(np.radians(self.vza_mean))

        # surfwater object:
        surfwater = xr.ones_like(self.raster.bands.isel(wl=0, drop=True).squeeze().astype(np.int8))
        surfwater.name = 'surfwater'
        self.raster['surfwater'] = surfwater

        # TODO remove or harmonize with landsat
        # self.U = self.raster.attrs['REFLECTANCE_CONVERSION_U']
        # convert into mW cm-2 um-1
        #self.solar_irradiance = xr.DataArray(self.raster.solar_irradiance / 10,
        #                                     coords={'wl': self.wl},
        #                                     attrs={
        #                                         'description': 'extraterrestrial solar irradiance from satellite metadata',
        #                                         'units': 'mW cm-2 um-1'})

        self.auxdatabase = auxdatabase
        self.output = output

        #########################
        # settings:
        #########################
        self.wl_process = self.raster.wl_to_process
        self.chunk = 512
        self._type = np.float32

        # self.sunglint_bands = [12]
        self.pressure_ref = 101500.
        self.iwl_swir = [-2, -1]
        self.bvis = 490
        self.bnir = 842
        self.bswir = 1610
        self.bswir2 = 2190

        # set extra bands (e.g., cirrus, water vapor)
        self.bcirrus = 1375
        self.cirrus = None
        self.bwv = 945
        self.wv = None

        # mask thresholding parameters
        self.sunglint_threshold = 0.2
        self.ndwi_threshold = 0.0
        self.vis_swir_index_threshold = 0.
        self.hcld_threshold = 3e-3

        # pre-computed auxiliary data
        self.dirdata = config['path']['grsdata']
        self.abs_gas_file = files('grs.data.lut.gases').joinpath('lut_abs_opt_thickness_normalized.nc')
        # self.lut_file = opj(self.dirdata, 'lut', 'opac_osoaa_lut_v2.nc')
        self.water_vapor_transmittance_file = files('grs.data.lut.gases').joinpath('water_vapor_transmittance.nc')
        self.load_auxiliary_data()

        # set retrieved parameter unit (Rrs or Lwn); is passed to fortran module
        self.rrs = False
        if self.output == 'Rrs':
            self.rrs = True

        #########################
        # variables:
        #########################

        self.headerfile = ''
        self.l2_product = None
        self.aux = None

        self.name = ''
        self.description = ''

        self.nodata = -999.9

        self.sza = []
        self.sazi = []
        self.vza = []
        self.vazi = []
        self.wkt = []  # geographical extent in POLYGON format
        self.outfile = ''
        self.aeronetfile = 'no'
        self.pressure_ref = 1015.20
        self.pressure = 1015.2
        self.ssp = 1015.2
        self.aot550 = 0.1
        self.aot = []
        self.angstrom = 1
        self.fcoef = 0.5
        self.rot = []
        self.oot = []

        # parameters for MAJA masks (CLM-cloud and MG2-geophysical)
        # masks order must follow bit order
        # CLM masks
        self.clm_masks = ['CLM_cloud_shadow', 'CLM_opaque_cloud', 'CLM_cloud_from_blue',
                          'CLM_cloud_multitemp', 'CLM_thin_cloud', 'CLM_shadow',
                          'CLM_shadow_from_outer_cloud', 'CLM_cirrus']
        # MG2 masks
        self.mg2_masks = ['MG2_Water_Mask', 'MG2_Cloud_Mask_All_Cloud', 'MG2_Snow_Mask',
                          'MG2_Shadow_Mask', 'MG2_Topographical_Shadows_Mask', 'MG2_Hidden_Areas_Mask',
                          ]
        self.maja_masks = np.concatenate([self.clm_masks, self.mg2_masks])

    def load_auxiliary_data(self):
        '''
        Load look-up tables data for gas absorption and backgroud transmittance

        :return:
        '''

        # get LUT
        self.gas_lut = xr.open_dataset(self.abs_gas_file)
        # self.aero_lut = xr.open_dataset(self.lut_file)
        # convert wavelength in nanometer
        # self.aero_lut['wl'] = self.aero_lut['wl'] * 1000
        # self.aero_lut['wl'].attrs['description'] = 'wavelength of simulation (nanometer)'
        self.Twv_lut = xr.open_dataset(self.water_vapor_transmittance_file)

    def set_outfile(self, file):
        '''
        Set name (path) of outputimage

        :param file:
        :return:
        '''
        self.outfile = file

    def set_aeronetfile(self, file):
        '''
        Deprecated

        :param file:
        :return:
        '''
        self.aeronetfile = file

    def get_flag(self, product, flag_name):
        '''
        Deprecated: obsolete solution now handled by GRSdriver module

        Get binary flag raster `flag_name` from `product`

        :param product: ProductIO object
        :param flag_name: name of the flag to be loaded
        :return:
        '''
        # TODO implement bit masks
        w, h = self.width, self.height
        flag_raster = np.zeros((w, h), dtype=np.int32, order='F').T

        return flag_raster

    # TODO implement DEM fetching
    def get_elevation(self, source='Copernicus30m'):

        self.elevation = xr.DataArray(np.zeros((self.height, self.width)), name="dem", coords=dict(
            y=self.y,
            x=self.x),
                                      attrs=dict(
                                          description="Digital elevation model from " + source,
                                          units="m")
                                      )

    def mu_N(self,
             sza,
             vza,
             azi,
             monoview=False):
        '''
        Compute the normal angle to wave slopes that produce sunglint.
        Warning: azi: azimuth in rad for convention azi=180 when sun-sensenor in oppositioon

        :param sza: solar zenith angle in degree
        :param vza: viewing zenith angle in degree
        :param azi: relative azimuth bewteen sun and sensor
        :param monoview: Set sensor viewing  configuration:
            - monoview = True : same viewing angles for all the spectral bands
            - monoview = False : viewing angles depend on spectral band (e.g. Sentinel-2 images)

        :return: cosine of normal angle to sunglint wave facets
        '''
        vzar = np.radians(vza)
        azir = np.radians(azi)
        szar = np.radians(sza)
        cos_alpha = np.cos(azir) * np.sin(vzar) * np.sin(szar) + np.cos(vzar) * np.cos(szar)
        xmu_N = (np.cos(szar) + np.cos(vzar)) / np.sqrt(2 * (1 + cos_alpha))
        if monoview:
            xmu_N = xmu_N.transpose('y', 'x')
        else:
            # ensure similar shape as inputs
            xmu_N = xmu_N.transpose('wl', 'y', 'x')
        return xmu_N

    def p_slope(self,
                sza,
                vza,
                azi,
                sigma2=0.02,
                monoview=False):
        '''
        Compute propability of wave slopes producing sunglint.

        :param sza: solar zenith angle in degree
        :param vza: viewing zenith angle in degree
        :param azi: relative azimuth bewteen sun and sensor
        :param sigma2: mean square slope of the wave slope distribution
        :param monoview: Set sensor viewing  configuration:
            - monoview = True : same viewing angles for all the spectral bands
            - monoview = False : viewing angles depend on spectral band (e.g. Sentinel-2 images)

        :return:
        '''

        cosN = self.mu_N(sza, vza, azi,monoview=monoview)
        thetaN = np.arccos(cosN)
        # stats == 'cm_iso':
        # TODO check consitency between sigma2 and formulation
        # Pdist_ = 1. / (np.pi *2.* sigma2) * np.exp(-1./2 * np.tan(thetaN) ** 2 / sigma2)
        xp_slope = 1. / (np.pi * sigma2) * np.exp(- np.tan(thetaN) ** 2 / sigma2) / cosN ** 4
        if monoview:
            xp_slope = xp_slope.transpose('y', 'x')
        else:
            xp_slope = xp_slope.transpose('wl', 'y', 'x')
        return xp_slope
