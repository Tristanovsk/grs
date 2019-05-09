# coding=utf-8
import re
import numpy as np
from dateutil import parser

from esasnappy import GPF, jpy
from esasnappy import Product, ProductUtils, ProductIO, ProductData
from esasnappy import FlagCoding, String, Mask

from .config import *


class info:
    def __init__(self, product, sensordata, aerosol='cams_forecast', ancillary='cams_forecast', output='Rrs'):
        '''

        :param product:
        :param sensordata:
        :param aerosol:
        :param ancillary:
        :param output: set the unit of the retrievals:
                 * 'Lwn', normalized water-leaving radiance (in mW cm-2 sr-1 μm-1)
                 * 'Rrs', remote sensing reflectance (in sr-1)
                 {default: 'Rrs']
        '''
        self.processor = __package__ + ' ' + VERSION
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
        self.lutfine = os.path.join(lut_root,
                                    sensordata.lutname + 'osoaa_band_integrated_aot0.01_aero_' + aero + '_ws2_pressure1015.2.nc')
        aero = 'rg0.80_sig0.60'
        self.lutcoarse = os.path.join(lut_root,
                                      sensordata.lutname + 'osoaa_band_integrated_aot0.01_aero_' + aero + '_ws2_pressure1015.2.nc')

        # set path for CAMS/ECMWF dataset
        self.cams_folder = cams_folder

        # set retrieved parameter unit (Rrs or Lwn); is passed to fortran module
        self.rrs = False
        if self.output == 'Rrs':
            self.rrs = True

        #########################
        # variables:
        #########################
        self.product = product
        self.headerfile = ''
        self.l2_product = None
        self.aux = None
        self.date = ''
        self.width = 0
        self.height = 0
        self.name = ''
        self.description = ''
        self.band_names = ''
        self.N = 0
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

    def set_outfile(self, file):
        self.outfile = file

    def set_aeronetfile(self, file):
        self.aeronetfile = file

    def get_product_info(self):
        product = self.product
        self.width = product.getSceneRasterWidth()
        self.height = product.getSceneRasterHeight()
        self.name = product.getName()
        self.description = product.getDescription()
        # self.band_names = product.getBandNames()
        self.date = parser.parse(str(product.getStartTime()))

    def get_bands(self, band_names=['B1']):
        '''get wavelengths, bands, geometries'''
        product = self.product
        self.band_names = band_names
        self.N = len(band_names)
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

    def get_flag(self, flag_name, rownum):
        '''get binary flag raster of row `rownum`'''

        flag_raster = np.zeros(self.width, dtype=np.int32)
        flag = self.product.getMaskGroup().get(flag_name)
        flag = jpy.cast(flag, Mask)
        flag.readPixels(0, rownum, self.width, 1, flag_raster)
        return flag_raster

    def load_data(self, rownum):
        # --------------------------------
        # construct arrays
        # --------------------------------

        self.mask = np.zeros(self.width, dtype=np.uint8, order='F')
        self.flags = np.zeros(self.width, dtype=np.uint8, order='F')

        self.hcld = np.zeros(self.width, dtype=self.type, order='F')
        self.rs2 = np.ma.zeros((self.N, self.width), dtype=self.type)
        self.vza = np.ma.zeros((self.N, self.width), dtype=self.type)
        self.razi = np.ma.zeros((self.N, self.width), dtype=self.type)
        self.sza = np.ma.zeros(self.width, dtype=self.type)
        self.sazi = np.ma.zeros(self.width, dtype=self.type)
        self.muv = np.ma.zeros((self.N, self.width), dtype=self.type)

        # set NDWI band
        self.B[self.sensordata.NDWI_vis].readPixels(0, rownum, self.width, 1, self.rs2[self.sensordata.NDWI_vis])
        self.B[self.sensordata.NDWI_nir].readPixels(0, rownum, self.width, 1, self.rs2[self.sensordata.NDWI_nir])
        self.ndwi = np.array((self.rs2[self.sensordata.NDWI_vis] - self.rs2[self.sensordata.NDWI_nir]) / \
                             (self.rs2[self.sensordata.NDWI_vis] + self.rs2[self.sensordata.NDWI_nir]))
        self.ndwi_band.loadRasterData()
        self.ndwi_band.setPixels(0, rownum, self.width, 1, self.ndwi)

        # set ndwi mask
        ndwi_ = (self.ndwi < self.sensordata.NDWI_threshold[0]) | (self.ndwi > self.sensordata.NDWI_threshold[1])
        self.mask[ndwi_] = 2

        # TODO for now not needed since SNAP API load all data irrespective of valid expression
        # # set valid expression
        # validexp = ' & ndwi >'+str(self.sensordata.NDWI_threshold[0])
        # def addValidPixelExpression(b,validexp):
        #     former = b.getValidPixelExpression()
        #     if former != None:
        #         b.setValidPixelExpression(former+validexp)
        #     else:
        #         b.setValidPixelExpression(validexp)

        # --------------------------------
        # load data
        # --------------------------------
        # addValidPixelExpression(self.SZA,validexp)
        self.SZA.readPixels(0, rownum, self.width, 1, self.sza)
        # addValidPixelExpression(self.SAZI,validexp)
        self.SAZI.readPixels(0, rownum, self.width, 1, self.sazi)
        self.mu0 = np.cos(np.radians(self.sza))

        def mask_array(arr, mask):
            return np.ma.array(arr, mask=mask, fill_value=np.nan)

        for iband in range(self.N):
            # addValidPixelExpression(self.B[iband],validexp)
            self.B[iband].readPixels(0, rownum, self.width, 1, self.rs2[iband])

            # check for nodata pixels
            nodata = self.B[iband].getGeophysicalNoDataValue()
            nodata_ = (self.rs2[iband].data == nodata)

            self.rs2[iband] = np.ma.array(self.rs2[iband],
                                          mask=nodata_ | ndwi_, fill_value=np.nan)

            mask = self.rs2[iband].mask
            # pin nodata pixels
            self.mask[nodata_] = 1

            self.VZA[iband].readPixels(0, rownum, self.width, 1, self.vza[iband])
            self.vza[iband] = mask_array(self.vza[iband], mask)
            self.VAZI[iband].readPixels(0, rownum, self.width, 1, self.razi[iband])
            self.razi[iband] = mask_array(self.razi[iband], mask)

            # get relative azimuth in OSOAA convention (=0 when sat and sun in opposition)
            self.razi[iband] = (180. - (self.razi[iband] - self.sazi)) % 360
            # self.razi[iband] = np.array([j % 360 for j in self.razi[iband]])
            self.muv[iband] = np.cos(np.radians(self.vza[iband]))

            # convert (if needed) into TOA reflectance
            if 'LANDSAT' in self.sensor:
                self.rs2[iband] = self.rs2[iband] * np.pi / (self.mu0 * self.U * self.solar_irr[iband] * 10)

        # --------------------------------
        # set mask cloud mask and/or export L1 flags
        # --------------------------------
        # TODO export L1 flags waiting for snap bug to be solved (subset remove mask info,
        #  https://forum.step.esa.int/t/problems-with-selecting-masks-as-input-in-graph-builder/3494/7 )
        try:
            cloud = self.get_flag(self.sensordata.cloud_flag, rownum)
            cirrus = self.get_flag(self.sensordata.cirrus_flag, rownum)
            self.flags = self.flags + ((cloud << 6) + (cirrus << 7))
        except:
            pass

        try:
            if self.sensordata.shadow_flag != '':
                shadow = self.get_flag(self.sensordata.shadow_flag, rownum)
                self.flags = (self.flags + (shadow << 8))
        except:
            pass

        # set high cloud cirrus mask
        try:
            self.product.getBand(self.sensordata.cirrus[0]).readPixels(0, rownum, self.width, 1, self.hcld)
            # convert (if needed) into TOA reflectance
            if 'LANDSAT' in self.sensor:
                self.hcld = self.hcld * np.pi / (self.mu0 * self.U * 366.97)
            self.flags = (self.flags + ((self.hcld > self.sensordata.cirrus[1]) << 5))
        except:
            pass  # print('No cirrus band available, high cloud flag discarded')

        # convert into FORTRAN 2D arrays (here, np.array)
        self.rs2 = np.array(self.rs2, order='F').T
        self.razi = np.array(self.razi, order='F').T
        self.vza = np.array(self.vza, order='F').T

    #     print('multiproc')
    #     with closing(Pool(8)) as p:
    #         return(p.map(self.f, range(self.N)))
    #
    # def f(self,iband):
    #     self.B[iband].readPixels(0, 0, self.width, self.height, self.rs2[iband])

    def unload_data(self):

        # unload data
        self.product.getBand('B10').unloadRasterData()
        self.SZA.unloadRasterData()
        self.SAZI.unloadRasterData()

        for i in range(self.N):
            self.B[i].unloadRasterData()
            self.VZA[i].unloadRasterData()
            self.VAZI[i].unloadRasterData()

    def create_product(self):

        product = self.product
        ac_product = Product('L2grs', 'L2grs', self.width, self.height)
        # writer = ProductIO.getProductWriter('NetCDF4-CF')  #
        writer = ProductIO.getProductWriter('BEAM-DIMAP')
        self.outfile_ext = self.outfile + '.dim'
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
        # TODO pass int or float arguement in proper format instead converting into string
        att_ = att('aerosol_data', ProductData.TYPE_ASCII)
        att_.setDataElems(self.aerosol)
        meta.addAttribute(att_)
        att_ = att('cams_file', ProductData.TYPE_ASCII)
        att_.setDataElems(self.aeronetfile)
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
        att_ = att('pressure', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.pressure))
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

        # set masks
        flags = ac_product.addBand('flags', ProductData.TYPE_UINT8)
        flags.setDescription('Flags for aquatic color purposes')
        # vflags = ac_product.addBand('valid', ProductData.TYPE_UINT8)
        # vflags.setDescription('used to set valid data pixels')

        # Also for each flag a layer should be created
        Color = jpy.get_type('java.awt.Color')
        coding = FlagCoding('flags')
        f0 = coding.addFlag("nodata", 1, "nodata in input image ")
        f1 = coding.addFlag("negative", 2, "negative values in visible ")
        f2 = coding.addFlag("ndwi", 4, "based on ndwi vis nir TOA ")
        f3 = coding.addFlag("ndwi_corr", 8, "based on ndwi vis nir after atmosperic correction ")
        f4 = coding.addFlag("high_nir", 16, "high radiance in the nir band (e.g., cloud, snow) ")
        f5 = coding.addFlag("hicld", 32, "high cloud as observed from cirrus band ")
        f6 = coding.addFlag("L1_cloud", 64, "opaque cloud flag from L1 image ")
        f7 = coding.addFlag("L1_cirrus", 128, "cirrus cloud flag from L1 image ")
        f8 = coding.addFlag("L1_shadow", 256, "cloud-shadow flag from L1 image ")

        ac_product.getFlagCodingGroup().add(coding)
        flags.setSampleCoding(coding)

        ac_product.addMask('mask_' + f0.getName(), 'flags.' + f0.getName(),
                           f0.getDescription(), Color.BLACK, 0.1)
        ac_product.addMask('mask_' + f1.getName(), 'flags.' + f1.getName(),
                           f1.getDescription(), Color.YELLOW, 0.1)
        ac_product.addMask('mask_' + f2.getName(), 'flags.' + f2.getName(),
                           f2.getDescription(), Color.RED, 0.1)
        ac_product.addMask('mask_' + f3.getName(), 'flags.' + f3.getName(),
                           f3.getDescription(), Color.PINK, 0.1)
        ac_product.addMask('mask_' + f4.getName(), 'flags.' + f4.getName(),
                           f4.getDescription(), Color.MAGENTA, 0.1)
        ac_product.addMask('mask_' + f5.getName(), 'flags.' + f5.getName(),
                           f5.getDescription(), Color.GREEN, 0.1)
        ac_product.addMask('mask_' + f6.getName(), 'flags.' + f6.getName(),
                           f6.getDescription(), Color.GRAY, 0.1)
        ac_product.addMask('mask_' + f7.getName(), 'flags.' + f7.getName(),
                           f7.getDescription(), Color.GRAY, 0.1)
        ac_product.addMask('mask_' + f8.getName(), 'flags.' + f8.getName(),
                           f8.getDescription(), Color.GRAY, 0.1)

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
            acband.setValidPixelExpression('mask_nodata == 0 && mask_ndwi == 0')
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
            acband.setValidPixelExpression('mask_nodata == 0 && mask_ndwi == 0')
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
        acband.setValidPixelExpression('mask_nodata == 0 && mask_ndwi == 0')
        ac_product.getBand(bname).setDescription('Glint reflection factor (BRDF) ')  # + self.band_names[iband])

        # estimated aerosol optical thickness at 550 nm
        bname = 'aot550'  # + self.band_names[iband]
        acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
        # acband.setSpectralWavelength(self.wl[iband])
        # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        acband.setValidPixelExpression('mask_nodata == 0 && mask_ndwi == 0')
        ac_product.getBand(bname).setDescription('aerosol optical thickness at 550 nm ')  # + self.band_names[iband])

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

        ac_product.setAutoGrouping(self.output + ':' + self.output + '_g_')

        ac_product.writeHeader(String(self.outfile_ext))
        # next line needed since snap 'writeHeader' force the extension to be consistent with data type (e.g., .tif for GeoTIFF)
        os.rename(self.outfile_ext, self.outfile_ext + '.incomplete')

        self.l2_product = ac_product

    def checksum(self, info):
        # TODO improve checksum scheme
        with open(self.outfile + '.checksum', "w") as f:
            f.write(info)

    def finalize_product(self):
        # TODO improve checksum scheme
        '''remove checksum file
        remove extension ".incomplete" from output file name'''
        os.remove(self.outfile + '.checksum')
        name = self.outfile_ext + '.incomplete'
        final_name = os.path.splitext(name)[0]
        os.rename(name, final_name)

    def print_info(self):
        ''' print info, can be used to check if object is complete'''
        print('Product: %s, %d x %d pixels, %s' % (self.name, self.width, self.height, self.description))
        print('Bands:   %s' % (list(self.band_names)))
        for i in range(len(self.wl)):
            print('Band ' + str(i) + ", " + str(self.B[i].getSpectralWavelength()) + ' centered on ' + str(
                self.wl[i]) + 'nm loaded')


class utils:

    def generic_resampler(self, s2_product, resolution=20, method='Bilinear'):
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
                  flag='FlagMedianAnd', opt=True):

        '''Resampling operator dedicated to Sentinel2-msi characteristics (e.g., viewing angles)
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
    def s2_resampler(product, resolution=20, upmethod='Bilinear', downmethod='First',
                     flag='FlagMedianAnd', opt=True):

        '''Resampling operator dedicated to Sentinel2-msi characteristics (e.g., viewing angles)
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

    def get_subset(self, product, wkt):
        '''subset from wkt POLYGON '''
        SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
        WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

        grid = WKTReader().read(wkt)
        op = SubsetOp()
        op.setSourceProduct(product)
        op.setGeoRegion(grid)
        op.setCopyMetadata(True)
        return op.getTargetProduct()

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
        ########
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

        return wkt

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
