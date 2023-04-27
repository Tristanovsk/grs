
import os, sys

import numpy as np
import pandas
from scipy.interpolate import interp1d
import netCDF4 as nc
import xarray as xr
import dask

import matplotlib.pyplot as plt

import logging
from dateutil import parser
import calendar, datetime

from . import config as cfg
from .utils import utils as u
from .acutils import aerosol
opj = os.path.join

class cams_product:
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

    def __init__(self,prod,
                 dir='./',
                 type='cams-global-atmospheric-composition-forecasts',
                 wls = [400, 440, 500, 550, 645, 670, 800, 865, 1020, 1240, 1640, 2130]):

        self.date = prod.date
        date_str = prod.date.strftime('%Y-%m-%d')
        self.dir = dir
        self.file = date_str + '-' + type + '.nc'
        xmin, ymin, xmax, ymax = prod.raster.rio.bounds()

        # lazy loading
        cams = xr.open_dataset(opj(dir, self.file), decode_cf=True,
                               chunks={'time': 1, 'x': 500, 'y': 500})
        # slicing
        cams = cams.sel(time=self.date, method='nearest')
        cams = cams.sel(latitude=slice(prod.latmax + 1, prod.latmin - 1),
                        longitude=slice(prod.lonmin - 1, prod.lonmax + 1))
        # rename "time" variable to avoid conflicts
        cams = cams.rename({'time':'time_cams'})
        if cams.u10.shape[0] == 0 or cams.u10.shape[1] == 0:
            print('no cams data, enlarge subset')

        self.raster = cams.interp(longitude=np.linspace(prod.lonmin, prod.lonmax, 12),
                   latitude=np.linspace(prod.latmax, prod.latmin, 12),
                   kwargs={"fill_value": "extrapolate"})
        self.raster= self.raster.rename({'longitude': 'x', 'latitude': 'y'})
        Nx = len(self.raster.x)
        Ny = len(self.raster.y)
        x = np.linspace(xmin, xmax, Nx)
        y = np.linspace(ymax, ymin, Ny)
        self.raster['x'] = x
        self.raster['y'] = y

        self.raster.rio.write_crs(prod.raster.rio.crs, inplace=True)

        param_ssa, param_aod = [], []
        for wl in wls:
            wl_ = str(wl)
            param_aod.append('aod' + wl_)
            param_ssa.append('ssa' + wl_)
        cams_aod = self.raster[param_aod].to_array(dim='wavelength')
        wl_cams = cams_aod.wavelength.str.replace('aod', '').astype(float)
        self.cams_aod = cams_aod.assign_coords(wavelength=wl_cams)
        cams_ssa = self.raster[param_ssa].to_array(dim='wavelength')
        wl_cams = cams_ssa.wavelength.str.replace('ssa', '').astype(float)
        self.cams_ssa = cams_ssa.assign_coords(wavelength=wl_cams)
        del cams_aod, cams_ssa

        self.variables = list(cams.keys())

        self.ozoneFactor = 46670.81518  # i.e., 300/6.428E-3
        self.pressureFactor = 0.01
        self.waterFactor = 0.1  # to convert in g/cm^2

        self.h2o = 3  # in g/cm-2
        # self.pressure = 1015.2  # in hPa
        self.o3du = 300  # in DU
        self.no2 = 0

        self.pressure_msl = 1015.2  # in hPa
        self.tco3 = np.nan
        self.tcwv = np.nan
        self.tcno2 = np.nan
        self.t2m = np.nan

        self.aot = []
        self.aot550 = 0.05
        self.aot_wl = []

        return

    def plot_params(self,params = ['aod550', 'aod2130', 'ssa550',
                                   't2m', 'msl', 'sp','tcco', 'tchcho',
                                   'tc_oh', 'tc_ch4', 'tcno2', 'gtco3',
                                   'tc_c3h8', 'tcwv', 'u10', 'v10']):

        Nrows = (len(params) + 1) // 4
        fig, axs = plt.subplots(Nrows, 4, figsize=(4 * 4.2, Nrows * 3.5))
        axs = axs.ravel()
        for i, param in enumerate(params):
            fig = self.raster[param].plot.imshow(robust=True, ax=axs[i])
            fig.axes.set_title(param)
            fig.colorbar.set_label(self.raster[param].units)
            fig.axes.set(xticks=[], yticks=[])
            fig.axes.set_ylabel('')
            fig.axes.set_xlabel('')
        plt.tight_layout()
        return fig, axs

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
            logging.info('downloading CAMS files...' + str(target))
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



    def subset_xr(self, ds, lonmin, lonmax, latmin, latmax, lat='latitude', lon='longitude'):
        '''

        :param ds:
        :param lonmin:
        :param lonmax:
        :param latmin:
        :param latmax:
        :param lat:
        :param lon:
        :return:
        '''
        mask_lon = (ds[lon] >= lonmin) & (ds[lon] <= lonmax)
        mask_lat = (ds[lat] >= latmin) & (ds[lat] <= latmax)
        return ds.where(mask_lon & mask_lat, drop=True)

    def get_xr_cams_aerosol(self, cams_file, product,
                            wls=[469, 550, 670, 865, 1240],
                            ):
        '''
        CAMS aerosol data loading, subset and Interpolation
        :param cams_file: absolute path of the CAMS netcdf file
        :param product: l2grs object
        :param wls: desired wavelengths to extract from database
        :return:
        '''

        N = len(wls)
        date = parser.parse(str(product.getStartTime()))
        wkt, lonmin, lonmax, latmin, latmax = u().get_extent(product)
        w, h = product.getSceneRasterWidth(), product.getSceneRasterHeight()

        self.aot = np.zeros(N, dtype=np.float32)
        self.aot_std = np.zeros(N, dtype=np.float32)
        self.aot_wl = wls

        # load CAMS netcdf file
        cams_xr = xr.open_dataset(cams_file)
        cams_xr = self.subset_xr(cams_xr, lonmin, lonmax, latmin, latmax)
        cams_xr = cams_xr.interp(time=date)

        for i, wl in enumerate(wls):
            param = 'aod' + str(wl)
            self.aot[i] = cams_xr[param].mean().data
            self.aot_std[i] = cams_xr[param].std().data
        self.aot550 = self.aot[1]
        self.aot550_std = self.aot_std[1]

        logging.info(f'{h}, {w}, {cams_xr.aod550.coords}')
        r, c = cams_xr.aod550.data.shape

        if (r > 1) and (c > 1):
            cams_rast = cams_xr.interp(longitude=np.linspace(lonmin, lonmax, w),
                                       latitude=np.linspace(latmax, latmin, h),
                                       kwargs={"fill_value": "extrapolate"})
            self.aot550rast = np.array(cams_rast['aod550'].data)
        else:
            self.aot550rast = np.full((w, h), self.aot550, order='F').T

        return  # u().getReprojected(prod, crs)

    def get_xr_cams_cds_aerosol(self, cams_file: str, l2h: object,
                                lutf: object, lutc: object,
                                wls=[400, 440, 500, 550, 645, 670, 800, 865, 1020, 1240, 1640, 2130],
                                ):
        '''
        CAMS aerosol data loading, subset and Interpolation
        :param cams_file: absolute path of the CAMS netcdf file
        :param l2h: l2grs object
        :param wls: desired wavelengths to extract from database
        :param i550: index of 550 nm wavelength in wls
        :return:
        '''

        wlsat, product = l2h.wl, l2h.product
        N = len(wlsat)
        date = parser.parse(str(product.getStartTime()))
        day = date.strftime(date.strftime('%Y-%m-%d'))
        wkt, lonmin, lonmax, latmin, latmax = u().get_extent(product)

        w, h = product.getSceneRasterWidth(), product.getSceneRasterHeight()

        param_ssa, param_aod = [], []
        for wl in wls:
            wl_ = str(wl)
            param_aod.append('aod' + wl_)
            param_ssa.append('ssa' + wl_)
        params = param_ssa + param_aod

        # ---------------------------------
        # open/load desired parameters
        # ---------------------------------
        cams_xr = xr.open_dataset(cams_file)[params]
        cams_daily = cams_xr.sel(time=day)

        # ---------------------------------
        # subset
        # ---------------------------------
        # first check longitude for nomenclature:
        # if close to meridian greenwich: -180 < lon < 180
        #   0 < lon < 360, otherwise
        # first convert from snap to CAMS nomenclature (0<lon<360)
        lonmin,lonmax=lonmin%360,lonmax%360
        if (lonmin > 350) and (lonmax < 10):
            lonmin = lonmin - 360
            cams_daily = cams_daily.assign_coords(longitude=(((cams_daily.longitude + 180) % 360) - 180)).sortby('longitude')

        lonslats = (lonmin, lonmax, latmin, latmax)

        # cams_sub = subset_xr(cams_daily, lonmin-1, lonmax+1, latmin-1, latmax+1)
        logging.info(f'cams cds  {lonmin}, {lonmax}, {latmin}, {latmax}')
        # cams_sub = self.subset_xr(cams_daily, lonmin, lonmax, latmin, latmax)
        cams_sub = cams_daily.interp(longitude=np.linspace(lonmin, lonmax, 12),
                                     latitude=np.linspace(latmax, latmin, 12),
                                     kwargs={"fill_value": "extrapolate"})

        # ---------------------------------
        # interpolate through dates
        # ---------------------------------
        cams_grs = cams_sub.interp(time=date, kwargs={"fill_value": "extrapolate"})

        # ---------------------------------
        # reshape to get ssa and aod as f(lon,lat,wavelength)
        # ---------------------------------
        cams_ssa = cams_grs[param_ssa].to_array(dim='wavelength')
        wl_cams = cams_ssa.wavelength.str.replace('ssa', '').astype(float)
        cams_ssa = cams_ssa.assign_coords(wavelength=wl_cams)
        cams_aod = cams_grs[param_aod].to_array(dim='wavelength')
        wl_cams = cams_aod.wavelength.str.replace('aod', '').astype(float)
        cams_aod = cams_aod.assign_coords(wavelength=wl_cams)

        aero = aerosol()

        r, c = cams_grs.aod550.shape
        # aot_550 = np.zeros((r, c), dtype=l2h.type)
        fcoef = np.zeros((r, c))
        aot_ = np.zeros((N, r, c))
        # ssa_grs = np.zeros((r, c, N), dtype=l2h.type)

        for ir in range(r):
            for ic in range(c):
                aero.fit_spectral_aot(cams_aod.wavelength, cams_aod.values[..., ir, ic])
                aot_[:, ir, ic] = aero.get_spectral_aot(wlsat)

        ssa_grs = cams_ssa.interp(wavelength=wlsat, kwargs={"fill_value": "extrapolate"}).chunk({'wavelength': N})
        aot_grs = ssa_grs.copy(data=aot_).chunk({'wavelength': N})
        aot_sca_grs = (ssa_grs * aot_grs).chunk({'wavelength': N})
        aot_sca_550 = aot_sca_grs.interp(wavelength=550, method='cubic')

        # normalization of Cext to get spectral dependence of fine and coarse modes
        nCext_f = lutf.Cext / lutf.Cext550
        nCext_c = lutc.Cext / lutc.Cext550

        for ir in range(r):
            for ic in range(c):
                fcoef[ir, ic] = \
                    aero.fit_aero(nCext_f, nCext_c, aot_sca_grs[..., ir, ic] / aot_sca_550[ir, ic])[0]
        fcoef = aot_sca_550.copy(data=fcoef)

        # self.ssa_grs = u().raster_regrid(ssa_grs, lonslats, h, w).values
        self.aot_grs = u().raster_regrid(aot_grs, lonslats, h, w).values
        self.aot_sca_grs = u().raster_regrid(aot_sca_grs, lonslats, h, w).values
        self.aot_sca_550 = u().raster_regrid(aot_sca_550, lonslats, h, w).values
        self.fcoef = u().raster_regrid(fcoef, lonslats, h, w).values

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
            logging.info('Error: not appropriate dataset for ecmwf/cams download')
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
            logging.info('Error: not appropriate cams settings for download')
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
            logging.error('Could not read file:', filename)
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
        logging.info(data.aot550, data.aot)
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
        logging.info(aot_df.columns)
        for i, col in enumerate(aot_df.columns):
            if 'AOD_' in col:
                wvs.append(int(col.replace('AOD_', '').replace('nm', '')))
                inds.append(i)

        wvs = np.array(wvs)
        inds = np.array(inds)
        logging.info(wvs)

        # wv_diffs = np.abs(wvs - 550)
        # logging.info(wv_diffs)
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
