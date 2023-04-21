import os, sys, re, glob
import numpy as np
import xarray as xr
import richdem as rd
from dateutil import parser
import logging
from pkg_resources import resource_filename

from . import config as cfg, auxdata, cams

opj = os.path.join


class product():
    '''

    :param product:
    :param sensordata:
    :param aerosol:
    :param ancillary:
    :param output: set the unit of the retrievals:
             * 'Lwn', normalized water-leaving radiance (in mW cm-2 sr-1 \mum-1)
             * 'Rrs', remote sensing reflectance (in sr-1)
             {default: 'Rrs']
    '''

    def __init__(self, l1c_obj=None, auxdatabase='cams', output='Rrs'):

        self.processor = __package__ + ' ' + cfg.VERSION

        self.raster = l1c_obj
        self.sensor = l1c_obj.attrs['satellite']
        self.date = self.raster.attrs['acquisition_date']
        self.width = self.raster.x.__len__()
        self.height = self.raster.y.__len__()
        self.vza_mean = self.raster.vza.mean()
        self.sza_mean = self.raster.vza.mean()
        self.air_mass_mean = 1. / np.cos(np.radians(self.sza_mean)) + 1. / np.cos(np.radians(self.vza_mean))

        self.sensordata = auxdata.sensordata(self.sensor)
        self.cams = cams.cams()
        self.U = l1c_obj.attrs['REFLECTANCE_CONVERSION_U']
        # convert into mW cm-2 um-1
        self.solar_irradiance = l1c_obj.solar_irradiance / 10
        self.auxdatabase = auxdatabase
        self.output = output

        #########################
        # settings:
        #########################
        self.block_size = 512
        self.sunglint_bands = [12]
        # data type for pixel values
        self.type = np.float32

        self.b565 = 560
        self.b865 = 865
        self.b1600 = 1610
        self.b2200 = 2190

        # mask thresholding parameters
        self.sunglint_threshold = 0.11
        self.ndwi_threshold = 0.01
        self.green_swir_index_threshold = 0.1
        self.hcld_threshold = 3e-3

        # pre-computed auxiliary data
        self.dirdata = resource_filename(__package__, '../grsdata/')
        self.abs_gas_file = opj(self.dirdata, 'gases', 'lut_abs_opt_thickness_normalized.nc')
        self.lut_file = opj(self.dirdata, 'lut', 'opac_osoaa_lut_v2.nc')
        self.water_vapor_transmittance_file = opj(self.dirdata, 'gases', 'water_vapor_transmittance.nc')
        self.load_auxiliary_data()

        # LUT for atmosphere radiance
        aero = 'rg0.10_sig0.46'
        self.lutfine = os.path.join(cfg.lut_root,
                                    self.sensordata.lutname +
                                    'osoaa_band_integrated_aot0.01_aero_' +
                                    aero + '_ws2_pressure1015.2.nc')
        aero = 'rg0.80_sig0.60'
        self.lutcoarse = os.path.join(cfg.lut_root,
                                      self.sensordata.lutname +
                                      'osoaa_band_integrated_aot0.01_aero_' +
                                      aero + '_ws2_pressure1015.2.nc')

        # set path for CAMS/ECMWF dataset
        self.cams_folder = cfg.cams_folder

        # set retrieved parameter unit (Rrs or Lwn); is passed to fortran module
        self.rrs = False
        if self.output == 'Rrs':
            self.rrs = True

        #########################
        # variables:
        #########################
        self.product = product
        # set band to be used in the process
        self.band_names = self.sensordata.band_names
        self.N = len(self.sensordata.band_names)

        self.headerfile = ''
        self.l2_product = None
        self.aux = None

        self.name = ''
        self.description = ''

        self.nodata = -999.9
        self.wl = []
        self.b = []
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

        # get LUT
        self.gas_lut = xr.open_dataset(self.abs_gas_file)
        # self.aero_lut = xr.open_dataset(self.lut_file)
        # convert wavelength in nanometer
        # self.aero_lut['wl'] = self.aero_lut['wl'] * 1000
        # self.aero_lut['wl'].attrs['description'] = 'wavelength of simulation (nanometer)'
        self.Twv_lut = xr.open_dataset(self.water_vapor_transmittance_file)

    def set_outfile(self, file):
        '''

        :param file:
        :return:
        '''
        self.outfile = file

    def set_aeronetfile(self, file):
        '''

        :param file:
        :return:
        '''
        self.aeronetfile = file

    def get_product_info(self):
        '''

        :return:
        '''
        product = self.product
        self.width = product.getSceneRasterWidth()
        self.height = product.getSceneRasterHeight()
        self.name = product.getName()
        self.description = product.getDescription()
        # self.band_names = product.getBandNames()
        print(product.getStartTime())
        if (((product.getStartTime())) is not None):
            self.date = parser.parse(str(product.getStartTime()))

    def get_flag(self, product, flag_name):
        '''
        get binary flag raster `flag_name` from `product`
        :param product: ProductIO object
        :param flag_name: name of the flag to be loaded
        :return:
        '''

        w, h = self.width, self.height
        flag_raster = np.zeros((w, h), dtype=np.int32, order='F').T
        flag = product.getMaskGroup().get(flag_name)
        flag = jpy.cast(flag, Mask)
        flag.readPixels(0, 0, w, h, flag_raster)
        return flag_raster

    def get_elevation(self, high_latitude=False):
        '''load elevation data into numpy array
        :param high_latitude: for |lat| > 60 deg, SRTM is not defined,
                if True, dem is set to GETASSE30
                '''

        self.elevation = np.zeros((self.width, self.height), dtype=self.type, order='F').T
        self.product = utils().add_elevation(self.product, high_latitude)
        dem = self.product.getBand('elevation')

        dem.readPixels(0, 0, self.width, self.height, self.elevation)
        self.elevation = np.array(self.elevation)

        return

    def load_data(self):
        '''
        load ta from input (subset) satellite image
        :return:
        '''
        # --------------------------------
        # construct arrays
        # --------------------------------
        w, h = self.width, self.height
        # set watermask as water (i.e., 1) for all pixels
        self.watermask = np.full((w, h), 1, dtype=np.uint8, order='F').T
        self.mask, self.flags = utils.init_fortran_array(2, (w, h), dtype=np.uint8)
        self.sza, self.sazi, arr = utils.init_fortran_array(3, (w, h))
        self.band_rad, self.vza, self.razi, self.muv = utils.init_arrayofarrays(4, [arr] * self.N)

        # --------------------------------
        # load Sun angles
        # --------------------------------
        self.SZA.readPixels(0, 0, w, h, self.sza)
        self.SAZI.readPixels(0, 0, w, h, self.sazi)
        self.mu0 = np.cos(np.radians(self.sza))

        # --------------------------------
        # load Masks
        # --------------------------------
        # get cirrus band if exists
        if self.sensordata.cirrus:
            logging.info(self.sensordata.cirrus[0])
            self.hcld = self.get_raster(self.product, self.sensordata.cirrus[0])
            # convert (if needed) into TOA reflectance
            if 'LANDSAT' in self.sensor:
                self.hcld = self.hcld * np.pi / (self.mu0 * self.U * 366.97)

        # get 02 band if exists
        if self.sensordata.O2band:
            self.O2band_raster = self.get_raster(self.product, self.sensordata.O2band[0])

        # if MAJA image provided, load (and write) AOT product band
        if self.maja:
            self.aot_maja = self.get_raster(self.maja, 'AOT_R1')
            self.l2_product.getBand('aot_maja').writePixels(0, 0, self.width, self.height, self.aot_maja)

        # --------------------------------
        # load data
        # --------------------------------
        # loop in reverse order to use the esasnappy "dispose()" function within the jvm
        for i, band in list(enumerate(self.band_names))[::-1]:
            logging.info(f'loading band {band}')
            self.B[i].readPixels(0, 0, w, h, arr)
            self.B[i].dispose()
            self.band_rad[i] = arr

            # check for nodata pixels and set mask
            nodata = self.B[i].getGeophysicalNoDataValue()
            nodata_ = self.band_rad[i] == nodata
            self.mask[nodata_] = 1

            if (self.mask == 1).all():
                raise ('No data available for the given subset of the image; process halted')

            # convert (if needed) into TOA reflectance
            if 'LANDSAT' in self.sensor:
                self.band_rad[i] = self.band_rad[i] * np.pi / (self.mu0 * self.U * self.solar_irr[i] * 10)

            self.VZA[i].readPixels(0, 0, w, h, arr)
            self.VZA[i].dispose()
            self.vza[i] = arr
            self.muv[i] = np.cos(np.radians(self.vza[i]))
            # get relative azimuth in OSOAA convention (=0 when sat and sun in opposition)
            self.VAZI[i].readPixels(0, 0, w, h, arr)
            self.VAZI[i].dispose()
            self.razi[i] = arr

            # convention RAZI = 0 when sun and satelite in opposition (Radiative transfer convention)
            self.razi[i] = (180. - (self.razi[i] - self.sazi)) % 360
            # convention RAZI = 180 when sun and satelite in opposition
        # self.razi[i] =  (self.razi[i] - self.sazi) % 360

        # self.razi[iband] = np.array([j % 360 for j in self.razi[iband]])

    def load_flags(self):
        '''
        Set flags from L1 data (must be called after `load_data`
        :return:
        '''

        # compute NDWI and set corresponding mask
        self.ndwi = np.array((self.band_rad[self.sensordata.NDWI_vis] - self.band_rad[self.sensordata.NDWI_nir]) /
                             (self.band_rad[self.sensordata.NDWI_vis] + self.band_rad[self.sensordata.NDWI_nir]))
        ndwi_ = (self.ndwi < self.sensordata.NDWI_threshold[0]) | (self.ndwi > self.sensordata.NDWI_threshold[1])
        self.mask[ndwi_] = 2
        self.flags = self.flags + (ndwi_ << 2)

        # compute NDWI and set corresponding mask
        self.ndwi_swir = np.array(
            (self.band_rad[self.sensordata.NDWI_swir_nir] - self.band_rad[self.sensordata.NDWI_swir_swir]) /
            (self.band_rad[self.sensordata.NDWI_swir_nir] + self.band_rad[self.sensordata.NDWI_swir_swir]))
        ndwi_swir_ = (self.ndwi_swir < self.sensordata.NDWI_swir_threshold[0]) | (
                self.ndwi_swir > self.sensordata.NDWI_swir_threshold[1])
        self.flags = self.flags + (ndwi_swir_ << 3)

        # --------------------------------
        # set mask cloud mask and/or export L1 flags
        # --------------------------------
        # TODO export L1 flags waiting for snap bug to be solved (subset remove mask info,
        #  https://forum.step.esa.int/t/problems-with-selecting-masks-as-input-in-graph-builder/3494/7 )
        # set high cloud cirrus mask
        if self.sensordata.cirrus:
            self.flags = self.flags + ((self.hcld > self.sensordata.cirrus[1]) << 5)
        if self.sensordata.O2band:
            self.flags = self.flags + ((self.O2band_raster > self.sensordata.O2band[1]) << 6)
        if self.sensordata.O2band:
            self.flags = self.flags + ((self.O2band_raster > self.sensordata.O2band[2]) << 7)

        try:
            cloud = self.get_flag(self.product, self.sensordata.cloud_flag)
            cirrus = self.get_flag(self.product, self.sensordata.cirrus_flag)
            self.flags = self.flags + (cloud << 8) + (cirrus << 9)
        except:
            pass

        if self.sensordata.shadow_flag != '':
            shadow = self.get_flag(self.product, self.sensordata.shadow_flag)
            self.flags = self.flags + (shadow << 10)

        # -------------------
        # for Sentinel 2
        # if MAJA L2A image is provided load MAJA flags
        mask_id = 11
        if self.maja:
            logging.info("copying MAJA masks...")
            # CLM masks
            masks = self.get_raster(self.maja, 'Aux_Mask_Cloud_R1', dtype=np.uint32)
            logging.info('CLM ' + str(np.unique(masks << mask_id)))
            logging.info('flags ' + str(np.unique(self.flags)))
            self.flags = self.flags + (masks << mask_id)
            logging.info('flags ' + str(np.unique(self.flags)))
            mask_id += len(self.clm_masks)

            # MG2 masks
            masks = self.get_raster(self.maja, 'Aux_Mask_MG2_R1', dtype=np.uint32)
            logging.info('MG2 ' + str(np.unique(masks << mask_id)))
            self.flags = self.flags + (masks << mask_id)
            logging.info('flags ' + str(np.unique(self.flags)))
            mask_id += len(self.mg2_masks)

        # -------------------
        # for Sentinel 2
        # if WaterDetect image is provided load Water Mask
        if self.waterdetect:
            water_true = (self.get_raster(self.waterdetect, 'band_1', dtype=np.uint32) == 1)
            self.flags = self.flags + (water_true << mask_id)
            # reinitialize watermask array to non-water
            self.watermask.fill(0)
            self.watermask[water_true] = 1

        return


def get_elevation(gdal_info_tgt, dem_glo30_dir, temp_dir=None, copy_dem_path=None):
    '''load elevation 10m data from Copernicus GLO30 into numpy array

            '''

    start_time_loc = datetime.utcnow()

    print('Computing DEM from GLO30 dir: %s' % dem_glo30_dir)

    # create temp dir

    if 'TMPDIR' in os.environ:

        temp_dir = os.environ['TMPDIR']

    else:

        temp_dir = os.path.abspath(os.getcwd())

    os.makedirs(temp_dir, exist_ok=True)

    temp_dir_session = None

    try:

        temp_dir_session = tempfile.mkdtemp(dir=temp_dir, prefix='grs2_dem_')

        # source band 10m from S2 L1C to read geometry. For example, B02 = 10m, B05=20m, B01=60m

        # try loading world vrt file containing info for all Copernicus GLO30 dems

        vrt_file = os.path.join(dem_glo30_dir, 'world.vrt')

        if not os.path.exists(vrt_file):
            # TODO : use gdalbuildvrt

            raise NotImplementedError('construction of .vrt file on the fly not implemented yet')

        # create .vrt file if general .vrt file does not exist

        target_file = os.path.join(temp_dir_session, 'dem.tif')

        print(' -> computing elevation from Copernicus GLO30 DEM, %s' % vrt_file)

        elevation = reproject_simple(vrt_file, gdal_info_tgt, target_file, resample_method='near', return_array=True)

        print(' -> min elevation=%s, max elevation=%s' % (np.min(elevation), np.max(elevation)))

        if copy_dem_path is not None:
            shutil.copy(target_file, copy_dem_path)

    finally:

        if temp_dir_session is not None:
            shutil.rmtree(temp_dir_session)

    print(' -> DEM generated in %s' % (datetime.utcnow() - start_time_loc))

    return elevation
