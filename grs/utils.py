# coding=utf-8
'''
Defines main python objects and image manipulation functions (linked to the ESA snappy library)
'''

import os, sys, re, glob
import numpy as np
import xarray as xr
import richdem as rd
from dateutil import parser
import logging

# from esasnappy import GPF, jpy
# from esasnappy import Product, ProductUtils, ProductIO, ProductData
# from esasnappy import FlagCoding, String, Mask

from . import config as cfg

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
        # TODO check if taking DEM files from cnes datalake is feasible -> faire un lien symbolique vers auxdata sur datalake
        op.setParameter("demName", "SRTM 3Sec")

        # srtm_path=cfg.srtm_path
        # logging.info(srtm_path)
        # for f in glob.glob(srtm_path+'/*.tif'):
        #    op.setParameter("externalDEMFile", f)

        if high_latitude:
            op.setParameter('demName', 'GETASSE30')
        op.setSourceProduct(product)
        return op.getTargetProduct()

    def get_dem_attributes(self, dem, sza=40, sun_azi=90):
        '''
        Compute slope, aspect and shades from dem for a given Sun geometry
        :param dem: Digital Elevation Model (numpy-like array)
        :param sza: Sun Zenith Angle in degrees
        :param sun_azi: Sun azimuth from North (clockwise)
        :return: slope in radians
        :return: shade for "hillshade values"
        '''
        sza = np.radians(sza)
        sun_azi = np.radians(sun_azi)
        dem = rd.rdarray(dem, no_data=np.nan)
        slope = np.radians(rd.TerrainAttribute(dem, attrib="slope_degrees"))
        aspect = rd.TerrainAttribute(dem, attrib="aspect")
        aspect_rd = np.radians(aspect)
        shade = np.cos(sza) * np.cos(slope) + np.sin(sza) \
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
            logging.error('sensor not recognized from input file')
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
            latitude=np.linspace(latmax, latmin, h), method="nearest",
            kwargs={"fill_value": "extrapolate"})

    @staticmethod
    def remove_na(arr):
        return arr[~np.isnan(arr)]
