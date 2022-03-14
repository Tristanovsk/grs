# coding=utf-8
'''
Defines main python objects and image manipulation functions (linked to the ESA snappy library)
'''

import os, sys, re, glob
import numpy as np
import xarray as xr
import richdem as rd
from dateutil import parser

from esasnappy import GPF, jpy
from esasnappy import Product, ProductUtils, ProductIO, ProductData
from esasnappy import FlagCoding, String, Mask

from . import config as cfg


class info:
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

    def __init__(self, product, sensordata, aerosol='cams_forecast', ancillary='cams_forecast', output='Rrs'):

        self.processor = __package__ + ' ' + cfg.VERSION
        self.sensordata = sensordata
        self.sensor = sensordata.sensor
        self.aerosol = aerosol
        self.ancillary = ancillary
        self.output = output

        #########################
        # settings:
        #########################
        # LUT for atmosphere radiance
        aero = 'rg0.10_sig0.46'
        self.lutfine = os.path.join(cfg.lut_root,
                                    sensordata.lutname + 'osoaa_band_integrated_aot0.01_aero_' + aero + '_ws2_pressure1015.2.nc')
        aero = 'rg0.80_sig0.60'
        self.lutcoarse = os.path.join(cfg.lut_root,
                                      sensordata.lutname + 'osoaa_band_integrated_aot0.01_aero_' + aero + '_ws2_pressure1015.2.nc')

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
        self.date = ''
        self.width = 0
        self.height = 0
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
        self.solar_irr = []
        self.U = 1.

        # data type for pixel values
        self.type = np.float32

        ### Define filtering thresholds
        self.hcld_threshold = 3e-3

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
        self.date = parser.parse(str(product.getStartTime()))

    def get_bands(self, band_names=None):
        '''
        get wavelengths, bands, geometries
        :param band_names:
        :return:
        '''

        product = self.product

        if not (band_names is None):
            self.band_names = band_names

        band_names = self.band_names
        self.N = len(self.band_names)
        N = range(self.N)
        self.VZA = [[] for i in N]
        self.VAZI = [[] for i in N]
        self.B = [[] for i in N]
        self.wl = [0 for i in N]
        # self.mvza = product.getBand('view_zenith_mean')
        # self.mazi = product.getBand('view_azimuth_mean')
        self.SZA = product.getBand('sun_zenith')
        self.SAZI = product.getBand('sun_azimuth')
        # self.SZA.loadRasterData()
        # self.SAZI.loadRasterData()

        for i in N:
            self.B[i] = product.getBand(band_names[i])
            # change for more accurate homemade values:
            # self.wl[i] = self.B[i].getSpectralWavelength()
            self.wl[i] = self.sensordata.central_wavelength[i]
            self.VZA[i] = product.getBand(self.sensordata.vza_name + band_names[i])
            self.VAZI[i] = product.getBand(self.sensordata.azi_name + band_names[i])
            # self.B[i].loadRasterData()
            # self.VZA[i].loadRasterData()
            # self.VAZI[i].loadRasterData()

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

    def get_raster(self, product, name, dtype=np.float):
        '''
        get geophysical raster `name` from `product`
        :param product: ProductIO object
        :param name: name of the band to be loaded
        :return:
        '''
        band = product.getBand(name)
        w, h = band.getRasterWidth(), band.getRasterHeight()
        raster = np.zeros((w, h), dtype=dtype, order='F').T
        band.readPixels(0, 0, w, h, raster)
        return raster

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
        # load Masks
        # --------------------------------
        # get cirrus band if exists
        if self.sensordata.cirrus:
            print(self.sensordata.cirrus[0])
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
        self.SZA.readPixels(0, 0, w, h, self.sza)
        self.SAZI.readPixels(0, 0, w, h, self.sazi)
        self.mu0 = np.cos(np.radians(self.sza))

        # loop in reverse order to use the esasnappy "dispose()" function within the jvm
        for i, band in list(enumerate(self.band_names))[::-1]:
            print('loading band ', band)
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
            # CLM masks
            masks = self.get_raster(self.maja, 'Aux_Mask_Cloud_R1', dtype=np.uint32)
            print('CLM', np.unique(masks << mask_id))
            print('flags', np.unique(self.flags))
            self.flags = self.flags + (masks << mask_id)
            print('flags', np.unique(self.flags))
            mask_id += len(self.clm_masks)

            # MG2 masks
            masks = self.get_raster(self.maja, 'Aux_Mask_MG2_R1', dtype=np.uint32)
            print('MG2', np.unique(masks << mask_id))
            self.flags = self.flags + (masks << mask_id)
            print('flags', np.unique(self.flags))
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

    def unload_data(self):
        '''unload data (not efficient due to ESA snappy issue with java jvm'''

        #
        self.product.getBand('B10').unloadRasterData()
        self.SZA.unloadRasterData()
        self.SAZI.unloadRasterData()

        for i in range(self.N):
            self.B[i].unloadRasterData()
            self.VZA[i].unloadRasterData()
            self.VAZI[i].unloadRasterData()

    def create_product(self, maja=None, waterdetect=None):
        '''
        Create output product dimensions, variables, attributes, flags....
        :return:
        '''

        # allocate flags param
        self.maja = maja
        self.waterdetect = waterdetect

        product = self.product
        ac_product = Product('L2grs', 'L2grs', self.width, self.height)
        writer = ProductIO.getProductWriter('NetCDF4-BEAM')  #
        # writer = ProductIO.getProductWriter('BEAM-DIMAP')
        self.outfile_ext = self.outfile + '.nc'  # dim'
        ac_product.setProductWriter(writer)

        ProductUtils.copyGeoCoding(product, ac_product)
        ProductUtils.copyMetadata(product, ac_product)
        ac_product.setStartTime(product.getStartTime())
        ac_product.setEndTime(product.getEndTime())

        # set metadata: ancillary data used for processing
        meta = jpy.get_type('org.esa.snap.core.datamodel.MetadataElement')
        att = jpy.get_type('org.esa.snap.core.datamodel.MetadataAttribute')
        # att(name=string,type=int), type: 41L->ascii; 12L->int32;

        meta = meta('L2')
        meta.setName('Ancillary Data')
        att_ = att('sensor', ProductData.TYPE_ASCII)
        att_.setDataElems(self.sensor)
        meta.addAttribute(att_)
        att_ = att('processor', ProductData.TYPE_ASCII)
        att_.setDataElems(self.processor)
        meta.addAttribute(att_)
        att_ = att('wkt', ProductData.TYPE_ASCII)
        att_.setDataElems(self.wkt)
        meta.addAttribute(att_)
        # TODO pass int or float argument in proper format instead converting into string
        att_ = att('aerosol_data', ProductData.TYPE_ASCII)
        att_.setDataElems(self.aerosol)
        meta.addAttribute(att_)
        att_ = att('cams_file', ProductData.TYPE_ASCII)
        att_.setDataElems(self.ancillary)
        meta.addAttribute(att_)
        att_ = att('aeronet_file', ProductData.TYPE_ASCII)
        att_.setDataElems(self.aeronetfile)
        meta.addAttribute(att_)
        att_ = att('aot_spectral', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aot))
        meta.addAttribute(att_)
        att_ = att('aot550', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aot550))
        meta.addAttribute(att_)
        att_ = att('pressure_sea_level', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.pressure_msl))
        att_.setUnit('hPa')
        meta.addAttribute(att_)
        att_ = att('ozone', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.o3du))
        att_.setUnit('DU')
        meta.addAttribute(att_)
        att_ = att('water_vapor', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.h2o))
        att_.setUnit('g cm-2')
        meta.addAttribute(att_)
        att_ = att('no2', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.tcno2))
        meta.addAttribute(att_)
        att_ = att('wind_module', ProductData.TYPE_ASCII)
        att_.setDataElems(str('...'))
        meta.addAttribute(att_)
        att_ = att('temperature_2m', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.t2m))
        meta.addAttribute(att_)

        ac_product.getMetadataRoot().addElement(meta)

        # -------------------------------------
        # set masks / flags
        # -------------------------------------

        flags = ac_product.addBand('flags', ProductData.TYPE_UINT32)
        flags.setDescription('Flags for aquatic color purposes')
        # vflags = ac_product.addBand('valid', ProductData.TYPE_UINT8)
        # vflags.setDescription('used to set valid data pixels')
        expr_valid_pixel = 'mask_nodata == 0 && mask_negative == 0'

        # Also for each flag a layer should be created
        Color = jpy.get_type('java.awt.Color')
        colors = [Color.BLUE, Color.YELLOW, Color.RED, Color.PINK, Color.MAGENTA, Color.GREEN, Color.GRAY] * 10
        coding = FlagCoding('flags')
        f=[None]*32
        mask_id = 0
        f[mask_id] = coding.addFlag("nodata", 2 ** mask_id, "nodata in input image ")
        mask_id += 1
        f[mask_id] = coding.addFlag("negative", 2 ** mask_id, "negative values in visible ")
        mask_id += 1
        f[mask_id] = coding.addFlag("ndwi", 2 ** mask_id, "based on ndwi vis nir TOA based on bands "+
                            self.band_names[self.sensordata.NDWI_vis]+" and "+
                            self.band_names[self.sensordata.NDWI_nir]+
                            " for range "+str(self.sensordata.NDWI_threshold))
        mask_id += 1
        f[mask_id] = coding.addFlag("ndwi_corr", 2 ** mask_id, "based on ndwi vis nir after atmosperic correction ")
        mask_id += 1
        f[mask_id] = coding.addFlag("high_nir", 2 ** mask_id, "high radiance in the nir band (e.g., cloud, snow); condition Rrs_g at " +
                            self.band_names[self.sensordata.high_nir[0]] + " greater than " + str(
                            self.sensordata.high_nir[1]))
        mask_id += 1
        f[mask_id] = coding.addFlag("hicld", 2 ** mask_id, "high cloud as observed from cirrus band; condition Rtoa at band " +
                            self.sensordata.cirrus[0] + " greater than " + str(self.sensordata.cirrus[1]))
        mask_id += 1
        f[mask_id] = coding.addFlag("moderate_cloud_risk_O2band", 2 ** mask_id,
                                    "moderate risk of bright cloud as observed from O2 band; condition Rtoa at band " +
                            self.sensordata.O2band[0] + " greater than " + str(self.sensordata.O2band[1]))
        mask_id += 1
        f[mask_id] = coding.addFlag("high_cloud_risk_O2band", 2 ** mask_id,
                                    "high risk of bright cloud as observed from O2 band; condition Rtoa at band " +
                            self.sensordata.O2band[0] + " greater than " + str(self.sensordata.O2band[2]))
        mask_id += 1
        f[mask_id] = coding.addFlag("L1_opaque_clouds", 2 ** mask_id, " flag from L1 image ")
        mask_id += 1
        f[mask_id] = coding.addFlag("L1_cirrus", 2 ** mask_id, "cirrus cloud flag from L1 image ")
        mask_id += 1
        f[mask_id] = coding.addFlag("L1_shadow", 2 ** mask_id, "cloud-shadow flag from L1 image ")

        for i_f in range(mask_id):
            ac_product.addMask('mask_' + f[i_f].getName(), 'flags.' + f[i_f].getName(),
                               f[i_f].getDescription(), colors[i_f], 0.3)


        # -------------------
        # for Sentinel 2
        # if MAJA L2A / WaterDetect provided, load respective flags
        # WARNING: mask_id must remain smaller than 32 (binary coding)
        mask_id += 1
        additional_f = []
        if self.maja:

            for i, mask in enumerate(self.maja_masks):
                print('Mask binary ', mask_id, 2 ** mask_id)
                additional_f.append(coding.addFlag(mask, 2 ** mask_id, 'Mask ' + mask + ' imported from MAJA chain'))
                mask_id += 1

        if self.waterdetect:
            print('waterdetect mask ', mask_id, 2 ** mask_id)
            additional_f.append(
                coding.addFlag('WaterDetect', 2 ** mask_id, 'Water mask imported from WaterDetect processing'))

        ac_product.getFlagCodingGroup().add(coding)
        flags.setSampleCoding(coding)


        # -------------------
        # for Sentinel 2
        # if MAJA L2A / WaterDetect provided, load respective flags
        # WARNING: mask_id must remain smaller than 32 (binary coding)

        mask_id = 0
        if self.maja:

            for i, mask in enumerate(self.maja_masks):
                f = additional_f[mask_id]
                ac_product.addMask('' + f.getName(), 'flags.' + f.getName(),
                                   f.getDescription(), colors[mask_id], 0.3)
                mask_id += 1

        if self.waterdetect:
            f = additional_f[mask_id]
            ac_product.addMask('mask_' + f.getName(), 'flags.' + f.getName(),
                               f.getDescription(), colors[mask_id], 0.3)

        # set data
        # Water-leaving radiance + sunglint
        for iband in range(self.N):
            bname = self.output + '_g_' + self.band_names[iband]
            acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
            acband.setSpectralWavelength(self.wl[iband])
            acband.setSpectralBandwidth(self.B[iband].getSpectralBandwidth())
            acband.setModified(True)
            acband.setNoDataValue(np.nan)
            acband.setNoDataValueUsed(True)
            acband.setValidPixelExpression(expr_valid_pixel)
            if self.output == 'Lwn':
                ac_product.getBand(bname).setDescription(
                    'Water-leaving plus sunglint normalized radiance (Lwn + Lg) in mW cm-2 sr-1 μm-1 at ' +
                    self.band_names[iband])
            else:
                ac_product.getBand(bname).setDescription(
                    'Water-leaving plus sunglint remote sensing reflectance (Rrs + Lg/F0) in sr-1 at ' +
                    self.band_names[iband])

        # Water-leaving radiance
        for iband in range(self.N):
            bname = self.output + '_' + self.band_names[iband]
            acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
            acband.setSpectralWavelength(self.wl[iband])
            acband.setSpectralBandwidth(self.B[iband].getSpectralBandwidth())
            acband.setModified(True)
            acband.setNoDataValue(np.nan)
            acband.setNoDataValueUsed(True)
            acband.setValidPixelExpression(expr_valid_pixel)
            if self.output == 'Lwn':
                ac_product.getBand(bname).setDescription(
                    'Normalized water-leaving radiance in mW cm-2 sr-1 μm-1 at ' + self.band_names[iband])
            else:
                ac_product.getBand(bname).setDescription(
                    'Remote sensing reflectance in sr-1 at ' + self.band_names[iband])
        # Sunglint reflection factor
        # for iband in range(self.N):
        bname = 'BRDFg'  # + self.band_names[iband]
        acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
        # acband.setSpectralWavelength(self.wl[iband])
        # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        acband.setValidPixelExpression(expr_valid_pixel + ' && ' + bname + ' >= 0')
        ac_product.getBand(bname).setDescription('Glint reflection factor (BRDF) ')  # + self.band_names[iband])

        # estimated aerosol optical thickness at 550 nm
        bname = 'aot550'  # + self.band_names[iband]
        acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
        # acband.setSpectralWavelength(self.wl[iband])
        # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        acband.setValidPixelExpression(expr_valid_pixel + ' && ' + bname + ' >= 0')
        ac_product.getBand(bname).setDescription('aerosol optical thickness at 550 nm ')  # + self.band_names[iband])

        # if maja option is on, copy aerosol optical thickness product from MAJA L2A image
        if self.maja:
            bname = 'aot_maja'  # + self.band_names[iband]
            acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
            # acband.setSpectralWavelength(self.wl[iband])
            # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
            acband.setModified(True)
            acband.setNoDataValue(np.nan)
            acband.setNoDataValueUsed(True)
            acband.setValidPixelExpression(expr_valid_pixel + ' && ' + bname + ' >= 0')
            ac_product.getBand(bname).setDescription(
                'AOT product from MAJA processing (L2A)')  # + self.band_names[iband])

        # Viewing geometry
        acband = ac_product.addBand('SZA', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('SZA').setDescription('Solar zenith angle in deg.')

        acband = ac_product.addBand('VZA', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('VZA').setDescription('Mean viewing zenith angle in deg.')

        acband = ac_product.addBand('AZI', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('AZI').setDescription('Mean relative azimuth angle in deg.')

        acband = ac_product.addBand('elevation', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        # TODO update the name of DEM data set used in function of the input (e.g., SRTM, GETASSE)
        ac_product.getBand('elevation').setDescription('elevation from SRTM 3Sec in meters')

        acband = ac_product.addBand('slope', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('slope').setDescription('terrain slope computed from DEM')

        acband = ac_product.addBand('shade', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('shade').setDescription('terrain shades computed from DEM and sun geometry')

        ac_product.setAutoGrouping(self.output + ':' + self.output + '_g_')

        ac_product.writeHeader(String(self.outfile_ext))
        # next line needed since snap 'writeHeader' force the extension to be consistent with data type (e.g., .tif for GeoTIFF)
        os.rename(self.outfile_ext, self.outfile_ext + '.incomplete')

        self.l2_product = ac_product

    def checksum(self, info):
        '''
        Save info on the current processing stage
        :param info:
        :return:
        '''
        # TODO improve checksum scheme
        with open(self.outfile + '.checksum', "w") as f:
            f.write(info)

    def deleteCache(self):
        S2CacheUtils = jpy.get_type('org.esa.s2tbx.dataio.cache.S2CacheUtils')
        S2CacheUtils.deleteCache()

    def finalize_product(self):
        '''remove checksum file
        remove extension ".incomplete" from output file name
        convert into netcdf (compressed) from gpt and ncdump/nco tool'''
        # jpy.destroy_jvm()
        self.product.dispose()
        self.l2_product.dispose()

        os.remove(self.outfile + '.checksum')
        name = self.outfile_ext + '.incomplete'
        # final_name = os.path.splitext(name)[0]
        os.rename(name, self.outfile_ext)  # final_name)

    def print_info(self):
        ''' print info, can be used to check if object is complete'''
        print('Product: %s, %d x %d pixels, %s' % (self.name, self.width, self.height, self.description))
        print('Bands:   %s' % (list(self.band_names)))
        for i in range(len(self.wl)):
            print('Band ' + str(i) + ", " + str(self.B[i].getSpectralWavelength()) + ' centered on ' + str(
                self.wl[i]) + 'nm loaded')


class utils:
    ''' utils for arrays, images manipulations'''

    def __init__(self):

        pass

    @staticmethod
    def init_arrayofarrays(number_of_array, dim):
        '''
        Initialize Nd array of Nd dimensions (given by dim)

        example:
            arr1, arr2 = init_array(2,(2,3,10))
        gives
            arr1.shape --> (2,3,10)
            arr2.shape --> (2,3,10)

        :param number_of_array: int
        :param dim: list

        :return: `number_of_array` numpy arrays of shape `dim`
        '''

        for i in range(number_of_array):
            yield np.array(dim)

    @staticmethod
    def init_fortran_array(number_of_array, dim, dtype=np.float32):
        '''
        Initialize Nd array of Nd dimensions (given by dim) in FORTRAN memory order
        (i.e., compliant with JAVA order provided by snappy)

        example:
            arr1, arr2 = init_array(2,(2,3,10))
        gives
            arr1.shape --> (2,3,10)
            arr2.shape --> (2,3,10)

        :param number_of_array: int
        :param dim: list
        :param dtype: numpy type for arrays
        :return: `number_of_array` numpy arrays of shape `dim`
        '''

        for i in range(number_of_array):
            yield np.zeros(dim, dtype=dtype, order='F').T

    def generic_resampler(self, s2_product, resolution=20, method='Nearest'):
        '''method: Nearest, Bilinear'''
        GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

        HashMap = jpy.get_type('java.util.HashMap')
        BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')

        parameters = HashMap()
        parameters.put('targetResolution', resolution)
        parameters.put('upsampling', method)
        parameters.put('downsampling', 'First')
        parameters.put('flagDownsampling', 'FlagMedianAnd')
        parameters.put('resampleOnPyramidLevels', True)

        return GPF.createProduct('Resample', parameters, s2_product)

    @staticmethod
    def resampler(product, resolution=20, upmethod='Bilinear', downmethod='First',
                  flag='FlagMedianAnd', opt=False):

        '''
        Resampling operator dedicated to Sentinel2-msi characteristics (e.g., viewing angles)

        :param product: S2-msi product as provided by esasnappy ProductIO.readProduct()
        :param resolution: target resolution in meters
        :param upmethod: interpolation method ('Nearest', 'Bilinear', 'Bicubic')
        :param downmethod: aggregation method ('First', 'Min', 'Max', 'Mean', 'Median')
        :param flag: method for flags aggregation ('First', 'FlagAnd', 'FlagOr', 'FlagMedianAnd', 'FlagMedianOr')
        :param opt: resample on pyramid levels (True/False)
        :return: interpolated target product'''

        resampler = jpy.get_type('org.esa.snap.core.gpf.common.resample.ResamplingOp')

        op = resampler()
        op.setParameter('targetResolution', resolution)
        op.setParameter('upsampling', upmethod)
        op.setParameter('downsampling', downmethod)
        op.setParameter('flagDownsampling', flag)
        op.setParameter('resampleOnPyramidLevels', opt)
        op.setSourceProduct(product)

        return op.getTargetProduct()

    @staticmethod
    def s2_resampler(product, resolution=20, upmethod='Bilinear', downmethod='Mean',
                     flag='FlagMedianAnd', opt=False):

        '''
        Resampling operator dedicated to Sentinel2-msi characteristics (e.g., viewing angles)

        :param product: S2-msi product as provided by esasnappy ProductIO.readProduct()
        :param resolution: target resolution in meters (10, 20, 60)
        :param upmethod: interpolation method ('Nearest', 'Bilinear', 'Bicubic')
        :param downmethod: aggregation method ('First', 'Min', 'Max', 'Mean', 'Median')
        :param flag: method for flags aggregation ('First', 'FlagAnd', 'FlagOr', 'FlagMedianAnd', 'FlagMedianOr')
        :param opt: resample on pyramid levels (True/False)
        :return: interpolated target product'''

        res = str(resolution)

        resampler = jpy.get_type('org.esa.s2tbx.s2msi.resampler.S2ResamplingOp')

        op = resampler()
        op.setParameter('targetResolution', res)
        op.setParameter('upsampling', upmethod)
        op.setParameter('downsampling', downmethod)
        op.setParameter('flagDownsampling', flag)
        op.setParameter('resampleOnPyramidLevels', opt)
        op.setSourceProduct(product)

        return op.getTargetProduct()

    def add_elevation(self, product, high_latitude=False):
        '''

        :param product: product open with snappy
        :param high_latitude: for |lat| > 60 deg, SRTM is not defined,
                if True, dem is set to GETASSE30
        :return:
        '''

        addelevation = jpy.get_type('org.esa.snap.dem.gpf.AddElevationOp')

        op = addelevation()
        op.setParameterDefaultValues()
        # TODO check if taking DEM files from cnes datalake is feasible
        # op.setParameter("demName", "External DEM")
        #
        # srtm_path=cfg.srtm_path
        # print(srtm_path)
        # for f in glob.glob(srtm_path+'/*.tif'):
        #       op.setParameter("externalDEMFile", f)

        if high_latitude:
            op.setParameter('demName', 'GETASSE30')
        op.setSourceProduct(product)
        return op.getTargetProduct()

    def get_dem_attributes(self,dem,sza=40,sun_azi=90):
        '''
        Compute slope, aspect and shades from dem for a given Sun geometry
        :param dem: Digital Elevation Model (numpy-like array)
        :param sza: Sun Zenith Angle in degrees
        :param sun_azi: Sun azimuth from North (clockwise)
        :return: slope in radians
        :return: shade for "hillshade values"
        '''
        sza=np.radians(sza)
        sun_azi=np.radians(sun_azi)
        dem= rd.rdarray(dem,no_data=np.nan)
        slope= np.radians(rd.TerrainAttribute(dem, attrib="slope_degrees"))
        aspect = rd.TerrainAttribute(dem, attrib="aspect")
        aspect_rd = np.radians(aspect)
        shade = np.cos(sza) * np.cos(slope) + np.sin(sza)\
                 * np.sin(slope) * np.cos(sun_azi - aspect_rd)
        return slope, shade

    def get_subset_old(self, product, wkt):
        '''subset from wkt POLYGON '''
        SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
        WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

        grid = WKTReader().read(wkt)
        op = SubsetOp()
        op.setSourceProduct(product)
        op.setGeoRegion(grid)
        op.setCopyMetadata(True)
        return op.getTargetProduct()

    def get_subset(self, product, wkt):
        HashMap = jpy.get_type('java.util.HashMap')
        parameters = HashMap()
        # parameters.put('bandNames',product.getBandNames()[0])
        parameters.put('geoRegion', wkt)
        parameters.put('copyMetadata', True)
        return GPF.createProduct('Subset', parameters, product)

    def print_array(self, arr):
        np.set_printoptions(threshold=np.nan)
        print(arr)

    def getMinMax(self, current, minV, maxV):
        if current < minV:
            minV = current
        if current > maxV:
            maxV = current
        return [minV, maxV]

    def get_extent(self, product):
        '''Get corner coordinates of the ESA SNAP product(getextent)

        # int step - the step given in pixels'''

        step = 1
        lonmin = 999.99

        GeoPos = ProductUtils.createGeoBoundary(product, step)

        lonmax = -lonmin
        latmin = lonmin
        latmax = lonmax

        for element in GeoPos:
            try:
                lon = element.getLon()
                [lonmin, lonmax] = self.getMinMax(lon, lonmin, lonmax)
            except (NameError):
                pass
            try:
                # TODO: separate method to get min and max
                lat = element.getLat()
                [latmin, latmax] = self.getMinMax(lat, latmin, latmax)
            except (NameError):
                pass
        wkt = 'POLYGON((' + str(lonmax) + ' ' + str(latmax) + ',' + str(lonmax) + ' ' \
              + str(latmin) + ',' + str(lonmin) + ' ' + str(latmin) + ',' + str(lonmin) + ' ' \
              + str(latmax) + ',' + str(lonmax) + ' ' + str(latmax) + '))'

        return wkt, lonmin, lonmax, latmin, latmax

    def getReprojected(self, product, crs='EPSG:4326', method='Bilinear'):
        '''Reproject a esasnappy product on a given coordinate reference system (crs)'''

        HashMap = jpy.get_type('java.util.HashMap')
        GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

        parameters = HashMap()
        parameters.put('crs', crs)
        parameters.put('resampling', method)

        return GPF.createProduct('Reproject', parameters, product)

    def get_sensor(self, file):
        '''
        Get sensor type from file name
        :param file: file in standard naming
        :return: sensor type
        '''
        file = os.path.basename(file)
        if ('S2A' in file):
            sensor = 'S2A'
        elif ('S2B' in file):
            sensor = 'S2B'
        elif ('LC08' in file) | bool(re.search(r"L[C,O]8", file)):
            sensor = 'LANDSAT_8'
        elif ('LE07' in file) | ('LE7' in file):
            sensor = 'LANDSAT_7'
        elif ('LT05' in file) | ('LT5' in file):
            sensor = 'LANDSAT_5'
        # TODO add to log file
        else:
            print('sensor not recognized from input file')
            sys.exit(-1)
        return sensor

    @staticmethod
    def raster_regrid(raster: xr.DataArray, lonslats: list, h: int, w: int):
        '''
        regrid cams raster onto satellite image grid of dim h,w
        TODO for the moment basic 2-step bilinear interpolation,
        TODO can be improved by using appropriate regridder, see: xESMF
        :param raster:
        :param h: height of image grid
        :param w: width of image grid
        '''
        lonmin, lonmax, latmin, latmax = lonslats
        return raster.interp(longitude=np.linspace(lonmin, lonmax, 12),
                             latitude=np.linspace(latmax, latmin, 12),
                             kwargs={"fill_value": "extrapolate"}).interp(
                             longitude=np.linspace(lonmin, lonmax, 512),
                             latitude=np.linspace(latmax, latmin, 512),
                             kwargs={"fill_value": "extrapolate"}).interp(
                             longitude=np.linspace(lonmin, lonmax, w),
                             latitude=np.linspace(latmax, latmin, h),method="nearest",
                             kwargs={"fill_value": "extrapolate"})
    @staticmethod
    def remove_na(arr):
        return arr[~np.isnan(arr)]
