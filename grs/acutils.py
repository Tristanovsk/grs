'''
Atmospheric Correction utilities to manage LUT and atmosphere parameters (aerosols, gases)
'''

import os, sys
import numpy as np
import xarray as xr
import logging
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import curve_fit


class lut:
    '''Load LUT FROM RT COMPUTATION (OSOAA_h)
    '''

    def __init__(self, band):

        self.smac_bands = []
        self.N = len(band)
        N = self.N
        self.lut_generator = "OSOAA_h"
        self.wl = []
        self.Cext = []
        self.Cext550 = 0
        self.Csca = []
        self.Csca550 = 0
        self.vza = []
        self.sza = []
        self.azi = []
        self.aot = []
        self.refl = []

    def load_lut(self, lut_file, ind_wl, aot=[0.01, 0.05, 0.1, 0.3, 0.5, 0.8], vza_max=20, reflectance=True):
        '''load lut calculated from OSOAA code

            Arguments:
                * ``lut_file`` -- netcdf file where lut data are stored for a given aerosol model
                * ``ind_wl`` -- indices of the desired central wavelength in ``lut_file``
                * ``aot`` -- aerosol optical thickness at 550 nm for which lut are loaded
                * ``vza_max`` -- load all lut data for vza <= vza_max
                * ``reflectance`` -- if true: return data in reflectance unit, return normalized radiane otherwise

            Construct object with:
                * ``aot`` -- aerosol optical thickness at 550 nm for which lut are loaded
                * ``Cext`` -- aerosol extinction coefficient (spectral)
                * ``Cext550`` -- aerosol extinction coefficient at 550 nm
                * ``vza`` -- viewing zenith angle (in deg)
                * ``sza`` -- solar zenith angle (in deg)
                * ``azi`` -- relative azimuth between sun and sensor (in opposition when azi = 180)
                * ``wl`` -- central wavelength of the sensor bands
                * ``refl`` -- Top-of-atmosphere reflectance (or normalized radiance if reflectance == False);
                                xarray of dims: [wl, sza, azi, vza]
              '''

        self.aot = aot

        Naot = len(aot)
        ok = 0
        for iaot in range(Naot):

            file = lut_file.replace('aot0.01', 'aot' + str(aot[iaot]))
            lut = Dataset(file, mode='r')
            self.Cext = lut.variables['Cext'][ind_wl]
            self.Cext550 = lut.variables['Cext550'][0]
            self.vza = lut.variables['vza'][:]
            self.sza = lut.variables['sza'][:]
            # azimuth in raditive transfer convention (0 deg whe n Sun and sensor in opposition)
            self.azi = lut.variables['azi'][:]
            self.wl = lut.variables['wavelength'][ind_wl]
            # shrink vza range (unused for S2)
            ind_vza = self.vza <= vza_max
            self.vza = self.vza[ind_vza]
            # allocate lut array
            if ok == 0:
                ok = 1
                nrad = np.zeros((Naot, len(self.wl), len(self.sza), len(self.azi), len(self.vza)))

            # fill in lut array
            nrad[iaot, :, :, :, :] = lut.variables['Istokes'][ind_wl, :, :, ind_vza]

        if reflectance:
            # convert into reflectance

            for i in range(len(self.sza)):
                nrad[:, :, i, :, :] = nrad[:, :, i, :, :] / np.cos(np.radians(self.sza[i]))
        # print(nrad.shape)
        self.refl = self._toxr(nrad)
        # reformat for remote sensing azimuth convention
        self.refl = self.refl.drop_isel(azi=-1)
        new_azi = (180 - self.refl.azi.values) % 360
        self.refl['azi'] = new_azi
        self.refl = xr.concat([self.refl, self.refl.sel(azi=0).assign_coords({'azi': 360})], dim='azi').sortby('azi')

    def _toxr(self, arr):
        # arr = np.array(arr)

        return xr.DataArray(arr,
                            dims=('aot', 'wl', 'sza', 'azi', 'vza'),
                            coords={'aot': self.aot,
                                    'wl': self.wl,
                                    'sza': self.sza,
                                    'azi': self.azi,
                                    'vza': self.vza})

    def interp_n_slice(self, sza_: np.array, vza_: np.array, azi_: np.array):
        '''
        Linear Interpolation of the lut array on the given angles
        '''

        self.refl = self.refl.interp(azi=azi_).interp(vza=vza_).interp(sza=sza_)

    def interp_lut(self, points, values, x):
        '''expected x dims: [[sza1, azi1, vza1],[sza2, azi2, vza2]...]'''
        from scipy.interpolate import interpn

        interp = np.ma.masked_invalid(interpn(points, values, x, bounds_error=False))

        return interp

    def plot_lut(self, vza, azi, values):

        spl = RectBivariateSpline(azi, vza, values)

        azi_ = np.linspace(0, max(azi), 360)
        vza_ = np.linspace(0, max(vza), 150)
        values_ = spl(azi_, vza_, grid=True)

        r, theta = np.meshgrid(vza_, np.radians(azi_))
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        # ax.contourf(theta,r, values)
        quadmesh = ax.pcolormesh(theta, r, values_)
        ax.grid(True)
        fig.colorbar(quadmesh, ax=ax)


class aerosol:
    '''
    aerosol parameters and parameterizations
    '''

    def __init__(self):
        self.aot550 = 0.1
        self.wavelengths = []
        self.aot = []
        self.wl = []
        self.ang = []
        self.o3du = 300
        self.no2du = 300
        self.fcoef = 0.5
        self.popt = []

    def func(self, lnwl, a, b, c):
        '''function for spectral variation of AOT'''

        return (a + b * lnwl + c * lnwl ** 2)

    def fit_spectral_aot(self, wl, aot):
        '''call to get fitting results on AOT data'''
        lnwl = np.log(wl)
        self.popt, pcov = curve_fit(self.func, lnwl, np.log(aot))

    def get_spectral_aot(self, wl):
        '''set aot for a given set of wavelengths'''
        lnwl = np.log(wl)
        return np.exp(self.func(lnwl, *self.popt))

    def func_aero(self, Cext, fcoef):
        '''function to fit spectral behavior of bimodal aerosols
         onto aeronet optical thickness'''
        return fcoef * Cext[0] + (1 - fcoef) * Cext[1]

    def fit_aero(self, nCext_f, nCext_c, naot):
        '''Call to get fine mode coefficient based on fitting on AOT data.

        Arguments:
          * ``nCext_f`` -- Normalized extinction coefficient of the fine mode aerosols
          * ``nCext_c`` -- Normalized extinction coefficient of the coarse mode aerosols
          * ``naot``    -- Normalized spectral aerosol optical thickness

        Return values:
          The mixing ratio of the fine mode aerosol

        Notes:
              .
            '''
        self.fcoef, pcov = curve_fit(self.func_aero, [nCext_f, nCext_c], naot)
        return self.fcoef


class cams_params:
    def __init__(self, name, resol):
        self.name = name
        self.resol = resol


class gases():
    def __init__(self):
        # atmosphere auxiliary data
        # TODO get them from CAMS
        self.pressure = 1010
        self.to3c = 6.5e-3
        self.tno2c = 3e-6
        self.tch4c = 1e-2
        self.psl = 1013
        self.coef_abs_scat = 0.3


class gaseous_transmittance(gases):

    def __init__(self, prod, cams):

        gases.__init__(self)
        self.xmin, self.ymin, self.xmax, self.ymax = prod.raster.rio.bounds()
        self.prod = prod
        self.cams = cams
        self.gas_lut = prod.gas_lut
        self.Twv_lut = prod.Twv_lut
        self.SRF = self.prod.raster.SRF
        self.air_mass_mean = self.prod.air_mass_mean
        self.pressure = cams.raster.sp * 1e-2
        self.coef_abs_scat = 0.3
        self.Tg_tot_coarse = None
        self.cams_gases = {'ch4': cams_params('tc_ch4', 4),
                           'no2': cams_params('tcno2', 7),
                           'o3': cams_params('gtco3', 4),
                           'h2o': cams_params('tcwv', 1), }

    def Tgas_background(self):
        gl = self.gas_lut
        pressure = self.pressure.round(1)
        self.ot_air = (gl.co + self.coef_abs_scat * gl.co2 +
                       self.coef_abs_scat * gl.o2 +
                       self.coef_abs_scat * gl.o4) / 1000

        wl_ref = gl.wl
        SRF_hr = self.prod.raster.SRF.interp(wl_hr=wl_ref.values)
        vals = np.unique(pressure)
        vals = vals[~np.isnan(vals)]
        if len(vals) == 1:
            vals = np.concatenate([vals, 1.2 * vals])
        Tg_raster = []
        for val in vals:
            Tg = np.exp(- self.air_mass_mean * self.ot_air * val)
            Tg = Tg.rename({'wl': 'wl_hr'})

            Tg_int = []
            for label, srf in SRF_hr.groupby('wl'):
                srf = srf.dropna('wl_hr').squeeze()
                Tg_ = Tg.sel(wl_hr=srf.wl_hr)
                wl_integr = Tg_.wl_hr.values

                Tg_ = np.trapz(Tg_ * srf, wl_integr) / np.trapz(srf, wl_integr)
                Tg_int.append(Tg_)
            Tg_raster.append(xr.DataArray(Tg_int, name='Ttot', coords={'wl': SRF_hr.wl.values}
                                          ).assign_coords({'pressure': val}))
        Tg_raster = xr.concat(Tg_raster, dim='pressure')
        return Tg_raster.interp(pressure=pressure).drop_vars(['pressure'])

    def Tgas(self, gas_name):

        renorm = 1.
        if gas_name == 'h2o':
            renorm = 0.4

        cams_gas = self.cams_gases[gas_name].name
        resol = self.cams_gases[gas_name].resol

        lut_abs = self.gas_lut[gas_name]
        rounded = renorm * self.cams.raster[cams_gas].round(resol)
        wl_ref = self.gas_lut.wl
        SRF_hr = self.prod.raster.SRF.interp(wl_hr=wl_ref.values)
        vals = np.unique(rounded)
        vals = vals[~np.isnan(vals)]
        if len(vals) == 1:
            vals = np.concatenate([vals, 1.2 * vals])
        Tg_raster = []
        for val in vals:
            Tg = np.exp(- self.air_mass_mean * lut_abs * val)
            Tg = Tg.rename({'wl': 'wl_hr'})

            Tg_int = []
            for label, srf in SRF_hr.groupby('wl'):
                srf = srf.dropna('wl_hr').squeeze()
                Tg_ = Tg.sel(wl_hr=srf.wl_hr)
                wl_integr = Tg_.wl_hr.values

                Tg_ = np.trapz(Tg_ * srf, wl_integr) / np.trapz(srf, wl_integr)
                Tg_int.append(Tg_)
            Tg_raster.append(xr.DataArray(Tg_int, name='Ttot', coords={'wl': SRF_hr.wl.values}
                                          ).assign_coords({'tc': val}))
        Tg_raster = xr.concat(Tg_raster, dim='tc')
        return Tg_raster.interp(tc=rounded)

    def get_gaseous_optical_thickness(self):
        gas_lut = self.gas_lut

        ot_o3 = gas_lut.o3 * self.to3c
        ot_ch4 = gas_lut.ch4 * self.tch4c
        ot_no2 = gas_lut.no2 * self.tno2c
        ot_air = (gas_lut.co + self.coef_abs_scat * gas_lut.co2 +
                  self.coef_abs_scat * gas_lut.o2 +
                  self.coef_abs_scat * gas_lut.o4) * self.pressure / 1000
        self.abs_gas_opt_thick = ot_ch4 + ot_no2 + ot_o3 + ot_air

    def get_gaseous_transmittance(self):
        Tg_tot = self.Tgas('ch4') * self.Tgas('no2') * \
                 self.Tgas('o3') * self.Tgas('h2o') * self.Tgas_background()

        # Tg_other = Tg_other.rename({'longitude': 'x', 'latitude': 'y'})
        # Nx = len(Tg_other.x)
        # Ny = len(Tg_other.y)
        # x = np.linspace(self.xmin, self.xmax, Nx)
        # y = np.linspace(self.ymax, self.ymin, Ny)
        # Tg_other['x'] = x
        # Tg_other['y'] = y
        self.Tg_tot_coarse = Tg_tot
        # TODO remove interp for the whole object and proceed with loop on spectral bands to save memory
        return Tg_tot  # .interp(x=self.prod.raster.x, y=self.prod.raster.y)

    def correct_gaseous_transmittance(self):

        return

    def get_gaseous_transmittance_old(self):

        self.get_gaseous_optical_thickness()
        wl_ref = self.gas_lut.wl  # .values
        SRF_hr = self.SRF.interp(wl_hr=wl_ref.values)
        Tg = np.exp(- self.air_mass_mean * self.abs_gas_opt_thick)
        Tg = Tg.rename({'wl': 'wl_hr'})

        Tg_int = []
        for label, srf in SRF_hr.groupby('wl'):
            srf = srf.dropna('wl_hr').squeeze()
            Tg_ = Tg.sel(wl_hr=srf.wl_hr)
            wl_integr = Tg_.wl_hr.values

            Tg_ = np.trapz(Tg_ * srf, wl_integr) / np.trapz(srf, wl_integr)
            Tg_int.append(Tg_)

        self.Tg_other = xr.DataArray(Tg_int, name='Ttot', coords={'wl': SRF_hr.wl.values})

    def other_gas_correction(self, raster_name='masked_raster', variable='Rtoa_masked'):
        raster = self.__dict__[raster_name]
        attrs = raster[variable].attrs
        if attrs.__contains__('other_gas_correction'):
            if attrs['other_gas_correction']:
                print('raster ' + raster_name + '.' + variable + ' is already corrected for other gases transmittance')
                print('set attribute other_gas_correction to False to proceed anyway')
                return
        if self.Tg_other is None:
            self.get_gaseous_transmittance(self.air_mass_mean)
        raster[variable] = raster[variable] / self.Tg_other
        raster[variable].attrs['other_gas_correction'] = True

    def water_vapor_correction(self, raster_name='coarse_masked_raster', variable='Rtoa_masked'):
        raster = self.__dict__[raster_name]
        attrs = raster[variable].attrs
        if attrs.__contains__('water_vapor_correction'):
            if attrs['other_gas_correction']:
                print('raster ' + raster_name + '.' + variable + ' is already corrected for water vapor transmittance')
                print('set attribute other_gas_correction to False to proceed anyway')
                return

        if self.Twv_raster is None:
            print('xarray of water vapor transmittance is not set, please run get_wv_transmittance_raster(tcwv_raster)')
            return
        raster[variable] = raster[variable] / self.Twv_raster
        raster[variable].attrs['water_vapor_correction'] = True

    def get_wv_transmittance_raster(self, tcwv_raster):
        tcwv_vals = tcwv_raster.tcwv.round(1)
        tcwvs = np.unique(tcwv_vals)
        tcwvs = tcwvs[~np.isnan(tcwvs)]
        # TODO improve for air_mass raster
        Twvs = self.Twv_lut.Twv.interp(air_mass=self.air_mass_mean).interp(tcwv=tcwvs, method='linear').drop('air_mass')
        self.Twv_raster = Twvs.interp(tcwv=tcwv_vals, method='nearest')


class misc:
    '''
    Miscelaneous utilities
    '''

    @staticmethod
    def get_pressure(alt, psl):
        '''Compute the pressure for a given altitude
           alt : altitude in meters (float or np.array)
           psl : pressure at sea level in hPa
           palt : pressure at the given altitude in hPa'''

        palt = psl * (1. - 0.0065 * np.nan_to_num(alt) / 288.15) ** 5.255
        return palt
