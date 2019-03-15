import numpy as np
import pandas
from scipy.interpolate import interp1d
import netCDF4 as nc
import calendar, datetime

from esasnappy import ProductIO

from .config import *
from .utils import utils as u


class sensordata:
    '''Dictionnaries of the auxilliary data, functions to be applied for the sensors data to process.

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
        * ``smac_band`` -- band names for SMAC files format
        * ``rg`` -- Cox-Munk Fresnel reflection factor ratio (R(wl)/Rref(2190nm) [Harmel et al., 2018]
        * ``angle_generator`` -- function to compute angles from metadata
    '''

    def __init__(self, sensor):
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
                    'smac_dir': os.path.join(smac_root, 'Coef_S2A_CONT_'),
                    'smac_bands': ('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'),
                    'cirrus_band': 'B10',
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
                    'ndwi_conf': [2, 7, NDWI_threshold]},

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
                    #TODO update for S2B after update of SMAC by CESBIO
                    'smac_dir': os.path.join(smac_root, 'Coef_S2A_CONT_'),
                    'smac_bands': ('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'),
                    'cirrus_band': 'B10',
                    'lut_name': 'S2B/lut_',
                    'rot': [0.23745233, 0.15521662, 0.09104522, 0.04485828, 0.03557663,
                            0.02918357, 0.02351827, 0.01829234, 0.01554304, 0.00127606, 0.00037652],
                    'rglint': [1.27878376, 1.26024821, 1.24195311, 1.22253644, 1.21708874, 1.21269296,
                               1.20815455, 1.20243389, 1.19906708, 1.11793355, 1.],
                    'ndwi_conf': [2, 7, NDWI_threshold]},

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
                          'ang_proc': os.path.join(grs_root, 'landsat_angles/TM/landsat_angles'),
                          'smac_dir': os.path.join(smac_root, 'coef_LANDSAT4_'),
                          'smac_bands': ('b1_CONT', 'b2_CONT', 'b3_CONT', 'b4_CONT', 'b5_CONT', 'b7_CONT'),
                          'cirrus_band': 'no',
                          'lut_name': 'L4/lut_L4_',
                          'rot': [0.16389, 0.084813, 0.046737, 0.01785, 0.001091, 0.0003578],
                          'rglint': [1.28037, 1.25722, 1.24089, 1.21893, 1.12406, 1.0],
                          'ndwi_conf': [1, 3, NDWI_threshold]},

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
                          'ang_proc': os.path.join(grs_root, 'landsat_angles/TM/landsat_angles'),
                          'smac_dir': os.path.join(smac_root, 'coef_LANDSAT5_'),
                          'smac_bands': ('b1_CONT', 'b2_CONT', 'b3_CONT', 'b4_CONT', 'b5_CONT', 'b7_CONT'),
                          'cirrus_band': 'no',
                          'lut_name': 'L5/lut_L5_',
                          'rot': [0.163419, 0.085201, 0.046514, 0.017945, 0.001098, 0.0003576],
                          'rglint': [1.28046, 1.257556, 1.24096, 1.21924, 1.12462, 1.0],
                          'ndwi_conf': [1, 3, NDWI_threshold]},

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
                          'ang_proc': os.path.join(grs_root, 'landsat_angles/TM/landsat_angles'),
                          'smac_dir': os.path.join(smac_root, 'coef_LANDSAT7_'),
                          'smac_bands': ('b1_CONT', 'b2_CONT', 'b3_CONT', 'b4_CONT', 'b5_CONT', 'b7_CONT'),
                          'cirrus_band': 'no',
                          'lut_name': 'L7/lut_L7_',
                          'rot': [0.174702, 0.091113, 0.046100, 0.018231, 0.001169, 0.0003645],
                          'rglint': [1.27891, 1.25551, 1.236695, 1.215597, 1.12477, 1.0],
                          'ndwi_conf': [1, 3, NDWI_threshold]},

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
                          'ang_proc': os.path.join(grs_root, 'landsat_angles/OLI/l8_angles'),
                          'smac_dir': os.path.join(smac_root, 'Coef_LANDSAT8_'),
                          'smac_bands': ('440_1', '490_1', '560_1', 'PAN_1', '660_1', '860_1', '1630_1', '2250_1'),
                          'cirrus_band': 'cirrus',
                          'lut_name': 'L8/lut_L8_',
                          'rot': [0.22899, 0.16786, 0.08999, 0.077445, 0.04785, 0.015507, 0.00128, 0.000366],
                          'rglint': [1.28637, 1.271187, 1.249180, 1.243706, 1.231680, 1.206419, 1.125028, 1.000000],
                          'ndwi_conf': [2, 5, NDWI_threshold]}
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
        self.smac_dir = info['smac_dir']
        self.smac_bands = np.array(info['smac_bands'])[idx]
        self.angle_processor = info['ang_proc']
        self.angle_names = np.array(info['ang_names'])[idx]
        self.vza_name = info['vza_name']
        self.azi_name = info['azi_name']
        self.band_names = np.array(info['band_names'])[idx]
        self.cirrus = info['cirrus_band']
        self.solar_irr = info['solar_irr']
        self.NDWI_vis, self.NDWI_nir, self.NDWI_threshold = info['ndwi_conf']


class cams:
    def __init__(self):
        '''
        Unit Conversion
        PWC (Precipitable Water Content), Grib Unit [kg/m^2]
        MSL (Mean Sea Level pressure),    Grib Unit [Pa]
        OZO (Ozone),                      Grib Unit [kg/m^2]

        calculation for Ozone according to R. Richter (20/1/2016):
        ----------------------------------------------------------
        GRIB_UNIT = [kg/m^2]
        standard ozone column is 300 DU (Dobson Units),
        equals to an air column of 3 mm at STP (standard temperature (0 degree C) and pressure of 1013 mbar).

        Thus, molecular weight of O3 (M = 48): 2.24 g (equals to 22.4 liter at STP)

        300 DU = 3 mm  (equals to (0.3*48 / 2.24) [g/m^2])
         = 6.428 [g/m^2] = 6.428 E-3 [kg/m^2]

        Example:

        ozone (kg/m^2  = 300 /6.428E-3 DU )
        pressure (Pa = 1/100 hPa)
        water vapor (kg/m^2 = 10^3/10^4 = 0.1 g/cm^2)

        '''
        self.tile_dir = ''
        self.aux_dir = ''
        self.ecmwf_file = ''
        self._nrTiles = 0

        self.ozoneFactor = 46670.81518  # i.e., 300/6.428E-3
        self.pressureFactor = 0.01
        self.waterFactor = 0.1 # to convert in g/cm^2

        self.h2o = 3  # in g/cm-2
        #self.pressure = 1015.2  # in hPa
        self.o3du = 300  # in DU
        self.no2 = 0

        self.msl = 1015.2 # in hPa
        self.tco3 = np.nan
        self.tcwv = np.nan
        self.tcno2 = np.nan
        self.t2m =np.nan

        self.aot = []
        self.aot550 = 0.05
        self.aot_wl = []

    def get_tile_dir(self, file):
        '''file : absolute path of S2 directory
           set full path of tile'''
        import fnmatch
        import os
        filemask = 'S2A_*_L1C_*'
        GRANULE = 'GRANULE'
        granuleDir = os.path.join(file, GRANULE)
        filelist = sorted(os.listdir(granuleDir))
        for tile in filelist:
            if fnmatch.fnmatch(tile, filemask) == False:
                continue
            self._nrTiles += 1
        if self._nrTiles > 1:
            # TODO implement for multiple tiles'
            print('several tiles are present, processing for single tile only')
            return False
        self.tile_dir = os.path.join(granuleDir, filelist[0])

    def get_aux_dir(self):
        self.aux_dir = os.path.join(self.tile_dir, 'AUX_DATA')

    def get_ecmwf_file(self, product):
        meta = product.getMetadataRoot().getElement('Level-1C_User_Product')
        ecmwf = meta.getElement('Auxiliary_Data_Info').getAttribute('ECMWF_DATA_REF')
        self.ecmwf_file = str(ecmwf.getData())

    def get_ecmwf_data(self):
        '''

        '''
        from osgeo.gdal_array import BandReadAsArray
        from osgeo import gdal
        import glob

        try:
            if (self.ecmwf_file == ''):
                self.ecmwf_file = glob.glob(self.aux_dir + '/*ECMWF*')[0]
            dataSet = gdal.Open(self.ecmwf_file)
            self.msl = BandReadAsArray(dataSet.GetRasterBand(2)) / 100
            self.o3du = BandReadAsArray(dataSet.GetRasterBand(3)) * self.ozoneFactor
        except:
            print('ERROR: reading ECMWF file!')
        return

    def load_cams_data(self, target, date, grid='0.125/0.125',
                       param='125.210/137.128/151.128/165.128/166.128/167.128/206.128/207.210/213.210/214.210/215.210/216.210',
                       data_type='cams_reanalysis'):

        ''' generate aerosol data from cams of ECMWF
            subset on the image grid (POLYGON wkt)
            reproject on the image coordinate reference system (crs)'''

        day = str(int(date.strftime('%d')))
        month = int(date.strftime('%m'))
        year = int(date.strftime('%Y'))

        # ------------download data
        if not os.path.isfile(str(target)):
            print('downloading CAMS files...' + str(target))
            startDate = '%04d%02d%02d' % (year, month, 1)
            numberOfDays = calendar.monthrange(year, month)[1]
            lastDate = '%04d%02d%02d' % (year, month, numberOfDays)
            now = datetime.datetime.now() - datetime.timedelta(days=2)
            # restricted access to data from the future
            # set lastDate as today date
            # TODO warning this will download incomplete data file for near real time data
            if now < datetime.datetime.strptime(lastDate, '%Y%m%d'):
                lastDate = datetime.datetime.strftime(now, '%Y%m%d')

            requestDates = startDate + '/TO/' + lastDate
            self.download_erainterim(str(target), requestDates, param=param,
                                     grid=grid, data_type=data_type)
        return

    def get_cams_aerosol(self, target, date, wkt,
                         param=['aod469', 'aod550', 'aod670', 'aod865', 'aod1240']):
        '''
        Nearest neighbor in time
        :param target:
        :param date:
        :param wkt:
        :return:
        '''
        # load netcdf file and get date-time index for further processing
        ncfile = nc.Dataset(target)
        time = ncfile.variables['time']

        dates = nc.num2date(time[:], time.units, time.calendar)
        ncfile.close()
        idx_time = np.abs([d - date for d in dates]).argmin(0)

        # set band names from cams:
        band_names = [a + '_time' + str(idx_time) for a in param]
        self.aot_wl = [469, 550, 670, 865, 1240]

        N = band_names.__len__()
        self.aot = np.zeros(N, dtype=np.float32)
        self.aot_std = np.zeros(N, dtype=np.float32)

        # ------------read/subset/collocate data
        prod = ProductIO.readProduct(str(target))
        #print(wkt)
        prod = u().get_subset(prod, wkt)

        h = prod.getBand(band_names[0]).getRasterHeight()
        w = prod.getBand(band_names[0]).getRasterWidth()
        aot_rast = np.zeros([N,w, h], dtype=np.float32)

        for i in range(N):
            prod.getBand(band_names[i]).loadRasterData()
            prod.getBand(band_names[i]).readPixels(0, 0, w, h, aot_rast[i,...])
            # aot_rast = np.ma.array(aot_rast,mask=  np.logical_or(aot_rast<0, aot_rast >2), fill_value=np.nan)
            self.aot[i] = aot_rast[i].mean()
            self.aot_std[i] = aot_rast[i].std()
            prod.getBand(band_names[i]).unloadRasterData()
        self.aot550 = self.aot[1]
        self.aot550_std = self.aot_std[1]
        self.aot_rast = aot_rast


        return  # u().getReprojected(prod, crs)

    def get_cams_ancillary(self, target, date, wkt,
                           param=['msl', 'tco3', 'tcwv', 'tcno2', 't2m']):
        '''
        Nearest neighbor in time
        :param target:
        :param date:
        :param wkt:
        :return:
        '''

        # load netcdf file and get date-time index for further processing
        ncfile = nc.Dataset(target)
        time = ncfile.variables['time']

        dates = nc.num2date(time[:], time.units, time.calendar)
        ncfile.close()
        idx_time = np.abs([d - date for d in dates]).argmin(0)

        # set band names from cams:
        band_names = [a + '_time' + str(idx_time) for a in param]
        N = band_names.__len__()

        # ------------read/subset/collocate data
        prod = ProductIO.readProduct(str(target))
        prod = u().get_subset(prod, wkt)

        h = prod.getBand(band_names[0]).getRasterHeight()
        w = prod.getBand(band_names[0]).getRasterWidth()
        data_ = np.zeros([N, w, h], dtype=np.float32)

        for i in range(N):
            prod.getBand(band_names[i]).readPixels(0, 0, w, h, data_[i,...])
            prod.getBand(band_names[i]).unloadRasterData()

        # TODO keep spatial resolution and export data as band product
        self.msl = data_[0].mean()
        self.tco3 = data_[1].mean()
        self.tcwv = data_[2].mean()
        self.tcno2 = data_[3].mean()
        self.t2m = data_[4].mean()

        self.msl = self.msl * self.pressureFactor
        self.o3du = self.tco3 * self.ozoneFactor
        self.h2o = self.tcwv * self.waterFactor

        return  # u().getReprojected(prod, crs)

    def get_cams_aerosol_old(self, target, date, wkt, crs, data_type):
        ''' generate aerosol data from cams of ECMWF
            subset on the image grid (POLYGON wkt)
            reproject on the image coordinate reference system (crs)'''
        from .utils import utils as u
        import calendar

        day = str(int(date.strftime('%d')))
        month = int(date.strftime('%m'))
        year = int(date.strftime('%Y'))

        # band names from cams:
        band_names = ['aod469', 'aod550', 'aod670', 'aod865', 'aod1240']
        band_names = [a + '_time' + day for a in band_names]
        self.aot_wl = [469, 550, 670, 865, 1240]
        N = band_names.__len__()
        self.aot = np.zeros(N, dtype=np.float32)

        # ------------download data
        if not target.is_file():
            print('downloading CAMS files...')
            startDate = '%04d%02d%02d' % (year, month, 1)
            numberOfDays = calendar.monthrange(year, month)[1]

            lastDate = '%04d%02d%02d' % (year, month, numberOfDays)

            requestDates = startDate + '/TO/' + lastDate
            self.download_erainterim(str(target), requestDates, data_type=data_type)

        # ------------read/subset/collocate data
        prod = ProductIO.readProduct(str(target))
        print(wkt)
        prod = u().get_subset(prod, wkt)
        h = prod.getBand(band_names[0]).getRasterHeight()
        w = prod.getBand(band_names[0]).getRasterWidth()
        aot_rast = np.zeros([w, h], dtype=np.float32)

        for i in range(N):
            prod.getBand(band_names[i]).loadRasterData()
            prod.getBand(band_names[i]).readPixels(0, 0, w, h, aot_rast)
            # aot_rast = np.ma.array(aot_rast,mask=  np.logical_or(aot_rast<0, aot_rast >2), fill_value=np.nan)
            self.aot[i] = aot_rast.mean()
            prod.getBand(band_names[i]).unloadRasterData()
        self.aot550 = self.aot[1]

        return  # u().getReprojected(prod, crs)

    def download_erainterim(self, target, date, time='00:00:00', grid='0.125/0.125',
                            param='137.128/151.128/206.210/207.210/213.210/214.210/215.210/216.210',
                            area=None, data_type=''):
        ''' This function open a connexion through an existing CAMS/ERAIterim account and download the requested data.


        Arguments:

            * ``target`` -- path name for output file
            * ``time`` -- time of forecast or reanalysis
            * ``date`` -- acquisition date (ex: '2015-01-31')
            * ``grid`` -- resolution in degrees (ex: '0.125/0.125')
            * ``param`` -- ID numbers for the requsted data
            * ``area`` -- option to restrict area (for global dataset None or '90/-180/-90/180')


        Notes:
            * ``server`` -- connexion au server via API ECMWF
            * ``step`` -- forecast time hours
        '''

        from .ecmwfapi import ECMWFDataServer

        server = ECMWFDataServer()
        if area is None:
            area = '90/-180/-90/180'

        step = '0'
        if data_type == 'cams_forecast':
            class_ = 'mc'
            dataset = 'cams_nrealtime'
            time = '00:00:00'
            step = '0/6/12/18'
            type = 'fc'


        elif data_type == 'cams_reanalysis':
            class_ = 'mc'
            dataset = 'cams_reanalysis'
            time = '00:00:00/06:00:00/12:00:00/18:00:00'
            type = 'an'

        elif data_type == 'interim':
            class_ = 'ei'
            dataset = 'interim'
            type = 'an'
        else:
            print('Error: not appropriate dataset for ecmwf/cams download')
            sys.exit()

        try:
            server.retrieve({ \
                'class': class_, \
                'dataset': dataset, \
                'date': date, \
                'grid': grid, \
                'levtype': 'sfc', \
                'param': param, \
                'step': step, \
                'stream': 'oper', \
                'time': time, \
                'type': type, \
                'format': 'netcdf',
                'area': area,
                'target': target})
        except:
            print('Error: not appropriate cams settings for download')
            sys.exit()

        server = None
        return


class Aeronet:
    '''Contains functions for importing AERONET measurements.
    Modifified from:
       Copyright 2012 Robin Wilson and contributors listed in the CONTRIBUTORS file.'''
    import pandas

    @classmethod
    def import_aeronet_data(cls, data, filename, time):
        '''Imports data from an AERONET data file to a given ``aeronet`` object.

        Arguments:
          * ``data`` -- Aeronet object to store the paramaters of interest (e.g., wavelenghts, aot, ozone)
          * ``filename`` -- The filename of the AERONET file described above
          * ``time`` -- The date and time of the simulation you want to run, used to choose the AERONET data which is closest
            in time. Provide this as a string in almost any format, and Python will interpret it. For example, ``'12/03/2010 15:39'``. When dates are ambiguous, the parsing routine will favour DD/MM/YY rather than MM/DD/YY.

        Return value:
          The function will return ``s`` with the ``aero_profile`` and ``aot550`` fields filled in from the AERONET data.

        Notes: Beware, this function returns ``data`` from ``filename`` for the closest measurement time in the limit
        of plus or minus 2 days from ``time``.
        '''

        # Load in the data from the file
        # TODO read and get exact wavelengths of CIMEL sensor from the AERONET file
        try:
            # header = pandas.read_csv(filename, skiprows=6,header=None,nrows=1)
            # df = pandas.read_csv(filename, skiprows=7, na_values=['N/A',-999.0], names = header.values[0].astype('str'))
            df = pandas.read_csv(filename, skiprows=6, na_values=['N/A', -999.0])
        except:
            print('Could not read file:', filename)
            return False
        # Parse the dates/times properly and set them up as the index
        df['Date(dd-mm-yyyy)'] = df['Date(dd-mm-yyyy)'].apply(cls._to_iso_date)
        df['timestamp'] = df.apply(lambda s: pandas.to_datetime(s['Date(dd-mm-yyyy)'] + ' ' + s['Time(hh:mm:ss)']),
                                   axis=1)
        df.index = pandas.DatetimeIndex(df.timestamp)

        given_time = time  # dateutil.parser.parse(time, dayfirst=True)

        df['timediffs'] = np.abs(df.timestamp - given_time).astype('timedelta64[ns]')
        # keep data within + or - 2 days
        df = df[df.timediffs < pandas.Timedelta('2 day')]

        # Get the AOT data
        data.wavelengths, data.aot550, data.aot = cls._get_aot(df)
        print(data.aot550, data.aot)
        data.o3du, data.no2du = cls._get_gaseous_abs(df)

        return

    @classmethod
    def _get_model_columns(cls, df):
        refr_ind = []
        refi_ind = []
        wvs = []
        radii_ind = []
        radii = []

        for i, col in enumerate(df.columns):
            if 'REFR' in col:
                refr_ind.append(i)
            elif 'REFI' in col:
                refi_ind.append(i)
                wv = int(col.replace('REFI', '').replace('(', '').replace(')', ''))
                wvs.append(wv)
            else:
                try:
                    rad = float(col)
                except:
                    continue
                radii_ind.append(i)
                radii.append(rad)

        return refr_ind, refi_ind, wvs, radii_ind, radii

    @classmethod
    def _to_iso_date(cls, s):
        '''Converts the date which is, bizarrely, given as dd:mm:yyyy to the ISO standard
        of yyyy-mm-dd.'''
        spl = s.split(':')
        spl.reverse()

        return '-'.join(spl)

    @classmethod
    def _get_aot(cls, df):
        '''Gets the AOT data from the AERONET dataset, choosing the AOT at the closest time
        to the time requested, and choosing the AOT measurement at the wavelength closest
        to 550nm.'''
        inds = []
        aot = []
        for i, col in enumerate(df.columns):
            if 'AOD_' in col:
                inds.append(i)

        inds.append(len(df.columns) - 1)
        inds = np.array(inds)

        # Remove the columns for AOT wavelengths with no data
        aot_df = df.ix[:, inds]

        aot_df = aot_df.dropna(axis=1, how='all')
        aot_df = aot_df.dropna(axis=0, how='any')

        wvs = []
        inds = []
        print(aot_df.columns)
        for i, col in enumerate(aot_df.columns):
            if 'AOD_' in col:
                wvs.append(int(col.replace('AOD_', '').replace('nm', '')))
                inds.append(i)

        wvs = np.array(wvs)
        inds = np.array(inds)
        print(wvs)

        # wv_diffs = np.abs(wvs - 550)
        # print(wv_diffs)
        # aot_col_index = wv_diffs.argmin()
        #
        # if (wv_diffs[aot_col_index] > 70):
        #     warnings.warn('Using AOT measured more than 70nm away from 550nm as nothing closer available - could cause inaccurate results.')

        rowind = aot_df.timediffs.idxmin()
        # aot550 = aot_df.ix[rowind, aot_col_index]
        aot = aot_df.ix[rowind, range(len(wvs))]
        aot_interp = interp1d(wvs, aot)
        aot550 = aot_interp(550.0)
        return wvs, aot550, np.array(pandas.to_numeric(aot))

    @classmethod
    def _get_gaseous_abs(cls, df):
        'ozone in dobson'
        o3col = []
        no2col = []
        inds = []
        for i, col in enumerate(df.columns):
            if 'Ozone' in col:
                o3col = i
            elif 'NO2' in col:
                no2col = i

        inds = len(df.columns) - 1
        inds = np.array([o3col, no2col, inds])

        # Remove the columns for AOT wavelengths with no data
        df = df.ix[:, inds]
        df = df.dropna(axis=1, how='all')
        df = df.dropna(axis=0, how='any')

        rowind = df.timediffs.idxmin()

        o3du = df.ix[rowind, 0]
        no2du = df.ix[rowind, 1]
        return 1000 * o3du, 1000 * no2du
