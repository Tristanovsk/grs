'''
Module dedicated to handle CAMS data with link to Copernicus API.
'''


import os, sys

import numpy as np
import pandas
from scipy.interpolate import interp1d
import xarray as xr

import matplotlib.pyplot as plt

import logging
import calendar, datetime
import cdsapi

opj = os.path.join


class CamsProduct:
    '''
    Unit Conversion:

    - PWC (Precipitable Water Content), Grib Unit [kg/m^2]
    - MSL (Mean Sea Level pressure),    Grib Unit [Pa]
    - OZO (Ozone),                      Grib Unit [kg/m^2]

    Calculation for Ozone according to R. Richter (20/1/2016):
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

    def __init__(self, prod,
                 cams_file=None,
                 dir='./',
                 type='forecast',
                 suffix=''
                 ):
        '''

        :param prod: l1c product from modeul product.prod
        :param cams_file:
        :param dir:
        :param type:
        :param suffix:
        '''

        self.prod = prod
        self.date = prod.time
        self.date_str = str(prod.time.dt.strftime('%Y-%m-%d').values)
        self.type = type
        self.wls = [469, 550, 670, 865, 1240]

        self.lonmin, self.latmin, self.lonmax, self.latmax = prod.rio.transform_bounds(4326)
        self.area = [self.latmax + 1,
                     self.lonmin - 1,
                     self.latmin - 1,
                     self.lonmax + 1]

        self.variable = [
            '10m_u_component_of_wind', '10m_v_component_of_wind',
            '2m_temperature',
            'mean_sea_level_pressure', 'surface_pressure',
            'ammonium_aerosol_optical_depth_550nm', 'black_carbon_aerosol_optical_depth_550nm',
            'dust_aerosol_optical_depth_550nm',
            'nitrate_aerosol_optical_depth_550nm', 'organic_matter_aerosol_optical_depth_550nm',
            'sea_salt_aerosol_optical_depth_550nm',
            'secondary_organic_aerosol_optical_depth_550nm', 'sulphate_aerosol_optical_depth_550nm',
            'total_aerosol_optical_depth_1240nm',
            'total_aerosol_optical_depth_469nm',
            'total_aerosol_optical_depth_550nm',
            'total_aerosol_optical_depth_670nm',
            'total_aerosol_optical_depth_865nm',
            'total_column_carbon_monoxide',
            'total_column_methane',
            'total_column_nitrogen_dioxide',
            'total_column_ozone', 'total_column_water_vapour']

        if cams_file:
            self.file = os.path.basename(cams_file)
            self.dir = os.path.dirname(cams_file)
            self.filepath = cams_file
        else:
            self.dir = dir
            if not os.path.exists(dir):
                os.makedirs(dir)
            self.file = self.date_str + '-' + type + suffix+'.nc'
            self.filepath = opj(self.dir, self.file)

    def cams_download(self):
        '''
        Autodownload from CAMS api of the cams netcdf data for the region-of-interest of the input image.

        :return:
        '''

        c = cdsapi.Client()

        c.retrieve(
            'cams-global-atmospheric-composition-forecasts',
            {
                'date': self.date_str + '/' + self.date_str,
                'type': self.type,
                'format': 'netcdf'
                ,
                'variable': self.variable,
                'time': ['00:00', '12:00'],
                'leadtime_hour': ['0', '3', '6', '9'],
                'area': self.area,
            },
            self.filepath)

        return

    def load(self):
        '''
        Lazy loading and then resmapling of the CAMS data for the region and date of interest.

        :return:
        '''

        # set geographic extents
        xmin, ymin, xmax, ymax = self.prod.rio.bounds()
        lonmin, latmin, lonmax, latmax = self.lonmin, self.latmin, self.lonmax, self.latmax

        if not os.path.exists(self.filepath):
            self.cams_download()

        # lazy loading
        cams = xr.open_dataset(self.filepath, decode_cf=True,
                               chunks={'time': -1, 'x': 500, 'y': 500})
        cams = cams.sel(latitude=slice(latmax + 1, latmin - 1))
        # check if image is on Greenwich meridian and adapt longitude convention
        if cams.longitude.min()>=0:
            if lonmin <= 0 and lonmax >= 0:

                    cams = cams.assign_coords({"longitude": (((cams.longitude + 180) % 360) - 180)}).sortby('longitude')
            else:
                # set longitude between 0 and 360 deg
                lonmin, lonmax, = lonmin % 360, lonmax % 360

        # slicing
        cams = cams.sel(longitude=slice(lonmin - 1, lonmax + 1)).load()

        # rename "time" variable to avoid conflicts
        # cams = cams.rename({'time':'time_cams'})
        if cams.u10.shape[0] == 0 or cams.u10.shape[1] == 0:
            print('no cams data, enlarge subset')

        # temporal interpolation
        if (cams.time.values[0] < self.date.values) & (cams.time.values[-1] > self.date.values):
            # if date within cams time range proceed to temporal interpolation
            cams = cams.interp(time=self.date)
        else:
            # otherwise get the nearest date
            cams = cams.sel(time=self.date,method='nearest')

        # spatial interpolation
        self.raster = cams.interp(longitude=np.linspace(lonmin, lonmax, 12),
                                  latitude=np.linspace(latmax, latmin, 12),
                                  kwargs={"fill_value": "extrapolate"})
        self.raster = self.raster.rename({'longitude': 'x', 'latitude': 'y'})
        Nx = len(self.raster.x)
        Ny = len(self.raster.y)
        x = np.linspace(xmin, xmax, Nx)
        y = np.linspace(ymax, ymin, Ny)
        self.raster['x'] = x
        self.raster['y'] = y

        self.raster.rio.write_crs(self.prod.rio.crs, inplace=True)

        param_aod = []
        for wl in self.wls:
            wl_ = str(wl)
            param_aod.append('aod' + wl_)

        cams_aod = self.raster[param_aod].to_array(dim='wl')

        wl_cams = cams_aod.wl.str.replace('aod', '').astype(float)
        self.cams_aod = cams_aod.assign_coords(wl=wl_cams)

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

    def plot_params(self, params=['amaod550', 'bcaod550', 'duaod550', 'niaod550',
                                  'omaod550', 'ssaod550', 'soaod550', 'suaod550',
                                  'aod550',
                                  't2m', 'msl', 'sp',
                                  'tcco', 'tc_ch4', 'tcno2', 'gtco3',
                                  'tcwv', 'u10', 'v10'],
                    **kwargs):
        '''
        Function to plot the cams data extracted for date and region of interest.

        :param params: parameters to plot
        :param kwargs: kwargs for matplotlib plotting
        :return: fig, axs
        '''

        Nrows = (len(params) + 4) // 4
        fig, axs = plt.subplots(Nrows, 4, figsize=(4 * 4.2, Nrows * 3.5))
        axs = axs.ravel()
        [axi.set_axis_off() for axi in axs]
        for i, param in enumerate(params):
            fig = self.raster[param].plot.imshow(robust=True, ax=axs[i], **kwargs)
            fig.axes.set_title(param)
            fig.colorbar.set_label(self.raster[param].units)
            fig.axes.set(xticks=[], yticks=[])
            fig.axes.set_ylabel('')
            fig.axes.set_xlabel('')
        plt.tight_layout()
        return fig, axs
