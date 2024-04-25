'''
Atmospheric Correction utilities to manage LUT and atmosphere parameters (aerosols, gases)
'''

import os, sys
import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

from numba import njit, prange


@njit()
def _getnearpos(array, value):
    '''
    Get index of nearest neighbor in an array.

    :param array: array to search in
    :param value: value you want the nearest index
    :return:
    '''
    idx = (np.abs(array - value)).argmin()
    return idx


class Rasterization:
    '''
    Class to numba compile python code with heavy loops.
    '''

    def __init__(self,
                 monoview=False):
        '''
        Set compiled function for:
            - monoview = True : same viewing angles for all the spectral bands
            - monoview = False : viewing angles depend on spectral band (e.g. Sentinel-2 images)

        :param monoview: True or False
        '''
        self.monoview = monoview
        if monoview:
            self.interp_Rlut_rayleigh = self._interp_Rlut_rayleigh_mono
            self.interp_Rlut = self._interp_Rlut_mono
        else:
            self.interp_Rlut_rayleigh = self._interp_Rlut_rayleigh
            self.interp_Rlut = self._interp_Rlut

    @staticmethod
    @njit()
    def _interp_Rlut_rayleigh_mono(szas, _sza,
                                   vzas, _vza,
                                   azis, _azi,
                                   Nwl, Ny, Nx, lut):
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        mus = np.cos(np.radians(_sza))
        for _iy in range(Ny):
            for _ix in range(Nx):
                if np.isnan(_sza[_iy, _ix]):
                    continue
                isza = _getnearpos(szas, _sza[_iy, _ix])
                iazi = _getnearpos(azis, _azi[_iy, _ix])
                ivza = _getnearpos(vzas, _vza[_iy, _ix])
                for _iwl in range(Nwl):
                    arr_lut[_iwl, _iy, _ix] = lut[_iwl, isza, ivza, iazi] / mus[_iy, _ix]
        return arr_lut

    @staticmethod
    @njit()
    def _interp_Rlut_mono(szas, _sza,
                          vzas, _vza,
                          azis, _azi,
                          aot_refs, _aot_ref,
                          Nwl, Ny, Nx, lut):
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        mus = np.cos(np.radians(_sza))
        for _iy in range(Ny):
            for _ix in range(Nx):
                if np.isnan(_sza[_iy, _ix]):
                    continue
                isza = _getnearpos(szas, _sza[_iy, _ix])
                iazi = _getnearpos(azis, _azi[_iy, _ix])
                ivza = _getnearpos(vzas, _vza[_iy, _ix])
                iaot_ref = _getnearpos(aot_refs, _aot_ref[_iy, _ix])
                for _iwl in range(Nwl):
                    arr_lut[_iwl, _iy, _ix] = lut[iaot_ref, _iwl, isza, ivza, iazi] / mus[_iy, _ix]
        return arr_lut

    @staticmethod
    @njit()
    def _interp_Rlut_rayleigh(szas, _sza,
                              vzas, _vza,
                              azis, _azi,
                              Nwl, Ny, Nx, lut):
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        mus = np.cos(np.radians(_sza))
        for _iy in range(Ny):
            for _ix in range(Nx):
                if np.isnan(_sza[_iy, _ix]):
                    continue
                isza = _getnearpos(szas, _sza[_iy, _ix])

                for _iwl in range(Nwl):
                    iazi = _getnearpos(azis, _azi[_iwl, _iy, _ix])
                    ivza = _getnearpos(vzas, _vza[_iwl, _iy, _ix])
                    arr_lut[_iwl, _iy, _ix] = lut[_iwl, isza, ivza, iazi] / mus[_iy, _ix]
        return arr_lut

    @staticmethod
    @njit()
    def _interp_Rlut(szas, _sza,
                     vzas, _vza,
                     azis, _azi,
                     aot_refs, _aot_ref,
                     Nwl, Ny, Nx, lut):
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        mus = np.cos(np.radians(_sza))
        for _iy in range(Ny):
            for _ix in range(Nx):
                if np.isnan(_sza[_iy, _ix]):
                    continue
                isza = _getnearpos(szas, _sza[_iy, _ix])
                iaot_ref = _getnearpos(aot_refs, _aot_ref[_iy, _ix])
                for _iwl in range(Nwl):
                    iazi = _getnearpos(azis, _azi[_iwl, _iy, _ix])
                    ivza = _getnearpos(vzas, _vza[_iwl, _iy, _ix])
                    arr_lut[_iwl, _iy, _ix] = lut[iaot_ref, _iwl, isza, ivza, iazi] / mus[_iy, _ix]
        return arr_lut

    @staticmethod
    @njit()
    def _interp_Tlut(szas, _sza,
                     aot_refs, _aot_ref,
                     Nwl, Ny, Nx, lut):
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        for _iy in range(Ny):
            for _ix in range(Nx):
                if np.isnan(_sza[_iy, _ix]):
                    continue
                isza = _getnearpos(szas, _sza[_iy, _ix])
                iaot_ref = _getnearpos(aot_refs, _aot_ref[_iy, _ix])
                for _iwl in range(Nwl):
                    arr_lut[_iwl, _iy, _ix] = lut[iaot_ref, _iwl, isza]
        return arr_lut

    @staticmethod
    @njit()
    def _interp_aotlut(aot_refs, _aot_ref,
                       Nwl, Ny, Nx, lut):
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

        for _iy in range(Ny):
            for _ix in range(Nx):
                iaot_ref = _getnearpos(aot_refs, _aot_ref[_iy, _ix])
                for _iwl in range(Nwl):
                    arr_lut[_iwl, _iy, _ix] = lut[iaot_ref, _iwl]
        return arr_lut

    @staticmethod
    @njit()
    def _multiplicate(arrwl, raster, arresult):
        Nwl, Ny, Nx = arresult.shape
        for iwl in range(Nwl):
            arresult[iwl] = arrwl[iwl] * raster
        return arresult


class Aerosol:
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


class CamsParams:
    def __init__(self, name, resol):
        self.name = name
        self.resol = resol


class Gases():
    '''
     Intermediate class to set parameters for absorbing gases.
    '''

    def __init__(self):
        # atmosphere auxiliary data
        # TODO get them from CAMS
        self.pressure = 1010
        self.to3c = 6.5e-3
        self.tno2c = 3e-6
        self.tch4c = 1e-2
        self.psl = 1013
        self.coef_abs_scat = 0.3


class GaseousTransmittance(Gases):
    '''
    Class containing functions to compute rasters of the direct transmittance of the absorbing gases.
    '''

    def __init__(self, prod, cams):

        Gases.__init__(self)
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
        self.cams_gases = {'ch4': CamsParams('tc_ch4', 4),
                           'no2': CamsParams('tcno2', 7),
                           'o3': CamsParams('gtco3', 4),
                           'h2o': CamsParams('tcwv', 1), }

    def Tgas_background(self):
        '''
        Compute direct transmittance for background absorbing gases: :math:`CO,\ O_2,\ O_4`

        :return:
        '''
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
            for label, srf in SRF_hr.groupby('wl', squeeze=False):
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
        '''
        Compute hyperspectral transmittance for a given absorbing gas and
        convolve it with the spectral respnse functions of the sattelite sensor.

        :param gas_name: name of the absorbing gas, choose between:
            - 'h2o'
            - 'o3'
            - n2o'
        :return: Gaseous transmittance for satellite bands
        '''
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
            for label, srf in SRF_hr.groupby('wl', squeeze=False):
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
        '''
        Get gaseous optival thickness from total column integrated concentration.
        :return:
        '''

        gas_lut = self.gas_lut

        ot_o3 = gas_lut.o3 * self.to3c
        ot_ch4 = gas_lut.ch4 * self.tch4c
        ot_no2 = gas_lut.no2 * self.tno2c
        ot_air = (gas_lut.co + self.coef_abs_scat * gas_lut.co2 +
                  self.coef_abs_scat * gas_lut.o2 +
                  self.coef_abs_scat * gas_lut.o4) * self.pressure / 1000
        self.abs_gas_opt_thick = ot_ch4 + ot_no2 + ot_o3 + ot_air

    def get_gaseous_transmittance(self):
        '''
        Get the final total gaseous transmittance.
        :return:
        '''
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
        '''
        Obsolete function
        :return:
        '''
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
        '''
        Correct for transmittance of high altitude gases
        (i.e., no coupling with low atmosphere scattering)

        :param raster_name:
        :param variable:
        :return:
        '''
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
        '''
        Correct for water vapor transmittance.

        :param raster_name:
        :param variable:
        :return:
        '''
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
        '''
        Get transmittance raster for correction of the image.

        :param tcwv_raster:
        :return:
        '''
        tcwv_vals = tcwv_raster.tcwv.round(1)
        tcwvs = np.unique(tcwv_vals)
        tcwvs = tcwvs[~np.isnan(tcwvs)]
        # TODO improve for air_mass raster
        Twvs = self.Twv_lut.Twv.interp(air_mass=self.air_mass_mean).interp(tcwv=tcwvs, method='linear').drop('air_mass')
        self.Twv_raster = Twvs.interp(tcwv=tcwv_vals, method='nearest')


class Misc:
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

    @staticmethod
    def transmittance_dir(aot, air_mass, rot=0):
        return np.exp(-(rot + aot) * air_mass)

    @staticmethod
    def air_mass(sza, vza):
        return 1 / np.cos(np.radians(vza)) + 1 / np.cos(np.radians(sza))
