'''
Atmospheric Correction utilities to manage LUT and atmosphere parameters (aerosols, gases)
'''

import os, sys
import numpy as np
import xarray as xr

from matplotlib import pyplot as plt
from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import curve_fit

from . import config as cfg


class lut:
    '''Load LUT FROM RT COMPUTATION (OSOAA_h)
    '''

    def __init__(self, band):

        self.smac_bands = []
        self.N = len(band)
        N = self.N
        self.lut_genearator = "OSOAA_h"
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
                * ``azi`` -- relative azimuth between sun and sensor (in opposition when azi = 0)
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
            self.azi = lut.variables['azi'][:]
            self.wl = lut.variables['wavelength'][ind_wl]
            # shrink vza range (unused for S2)
            ind_vza = self.vza <= vza_max
            self.vza = self.vza[ind_vza]
            # allocate lut array
            if ok == 0:
                ok = 1
                nrad = np.ndarray((Naot, len(self.wl), len(self.sza), len(self.azi), len(self.vza)))

            # fill in lut array
            nrad[iaot, :, :, :, :] = lut.variables['Istokes'][ind_wl, :, :, ind_vza]

            if reflectance:
                # convert into reflectance

                for i in range(len(self.sza)):
                    nrad[:, :, i, :, :] = nrad[:, :, i, :, :] / np.cos(np.radians(self.sza[i]))

            self.refl = self._toxr(nrad)

    def _toxr(self, arr):
        arr = np.array(arr)

        return xr.DataArray(arr,
                            dims=('aot', 'wl', 'sza', 'azi', 'vza'),
                            coords={'aot': self.aot,
                                    'wl': self.wl,
                                    'sza': self.sza,
                                    'azi': self.azi,
                                    'vza': self.vza})

    def interp_n_slice(self,sza_:np.array,vza_:np.array,azi_:np.array):
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


class smac:
    ''' Gaseous absorption and transmission from pre-calculated 6S/SMAC data '''

    def __init__(self, smac_bands, smac_dir):

        self.smac_bands = smac_bands
        self.N = len(smac_bands)
        N = self.N

        ##################
        # FROM SMAC CESBIO
        self.smacdir = os.path.join(cfg.smac_root, smac_dir)
        self.uh2o = 3  # Water vapour (g/cm2)
        self.ah2o = [0] * N  # coef from SMAC computed from 6S th2o  = np.exp ( (ah2o) * ( (uh2o*m)  ** (nh2o) ) )
        self.nh2o = [0] * N  # coef from SMAC computed from 6S
        # O3
        self.uo3 = 330 * 0.001  # in cm.atm (= DU/1000)
        self.ao3 = [0] * N
        self.no3 = [0] * N
        # O2
        self.ao2 = [0] * N
        self.no2 = [0] * N
        self.po2 = [0] * N
        # CO2
        self.aco2 = [0] * N
        self.nco2 = [0] * N
        self.pco2 = [0] * N
        # NH4
        self.ach4 = [0] * N
        self.nch4 = [0] * N
        self.pch4 = [0] * N
        # NO2
        self.ano2 = [0] * N
        self.nno2 = [0] * N
        self.pno2 = [0] * N
        # CO
        self.aco = [0] * N
        self.nco = [0] * N
        self.pco = [0] * N
        ##################

        self.tg = []  # gaseous transmittance (up and down)

    def set_values(self, o3du=300, h2o=0, no2=0):
        '''set atmospheric concentration values:
        :param o3du: ozone in DU
        :param h2o: water vapor in g/cm^2
        :param no2: nitrous dioxide in ...'''

        self.uo3 = o3du / 1000  # conversion DU to cm.atm
        self.uh2o = h2o
        self.uno2 = no2
        return

    # TODO write a function to set standard values for all compounds
    def set_standard_values(self, peq):
        i = 0
        # gaseous transmissions (downward and upward paths)
        self.uo2 = (peq ** (self.po2[i]))
        self.uco2 = (peq ** (self.pco2[i]))
        self.uch4 = (peq ** (self.pch4[i]))
        self.uno2 = (peq ** (self.pno2[i]))
        self.uco = (peq ** (self.pco[i]))

        return

    def set_gas_param(self):
        '''Load gaseous absorption parameters as computed for SMAC from 6S RT code'''
        for i in range(self.N):
            try:
                f = open(self.smacdir + self.smac_bands[i] + '.dat', 'r')
                lines = f.readlines()
                f.close()
                # H20
                temp = lines[0].strip().split()
                self.ah2o[i] = float(temp[0])
                self.nh2o[i] = float(temp[1])
                # O3
                temp = lines[1].strip().split()
                self.ao3[i] = float(temp[0])
                self.no3[i] = float(temp[1])
                # O2
                temp = lines[2].strip().split()
                self.ao2[i] = float(temp[0])
                self.no2[i] = float(temp[1])
                self.po2[i] = float(temp[2])
                # CO2
                temp = lines[3].strip().split()
                self.aco2[i] = float(temp[0])
                self.nco2[i] = float(temp[1])
                self.pco2[i] = float(temp[2])
                # NH4
                temp = lines[4].strip().split()
                self.ach4[i] = float(temp[0])
                self.nch4[i] = float(temp[1])
                self.pch4[i] = float(temp[2])
                # NO2
                temp = lines[5].strip().split()
                self.ano2[i] = float(temp[0])
                self.nno2[i] = float(temp[1])
                self.pno2[i] = float(temp[2])
                # CO
                temp = lines[6].strip().split()
                self.aco[i] = float(temp[0])
                self.nco[i] = float(temp[1])
                self.pco[i] = float(temp[2])
            except:
                print('WARNING, NO SMAC FILES FOUND FOR BAND ' + self.smac_bands[i] + ' in ' + self.smacdir + '!')
                sys.exit(0)

    def compute_gas_trans(self, iband, pressure, mu0, muv):
        '''Compute gaseous transmittances (up and down) from SMAC parameters and ancillary data
           pressure : actual pressure in hPA
           mu0 : cosine of solar zenith angle
           muv : cosine of viewing zenith angle'''

        peq = pressure / 1013.25

        # air mass
        m = [1 / mu0 + 1 / muv]

        i = iband
        # gaseous transmissions (downward and upward paths)
        self.uo2 = (peq ** (self.po2[i]))
        self.uco2 = (peq ** (self.pco2[i]))
        self.uch4 = (peq ** (self.pch4[i]))
        self.uno2 = (peq ** (self.pno2[i]))
        self.uco = (peq ** (self.pco[i]))

        to3 = np.array([np.exp((self.ao3[i]) * ((self.uo3 * x) ** (self.no3[i]))) for x in m])
        th2o = np.array([np.exp((self.ah2o[i]) * ((self.uh2o * x) ** (self.nh2o[i]))) for x in m])
        to2 = np.array([np.exp((self.ao2[i]) * ((self.uo2 * x) ** (self.no2[i]))) for x in m])
        tco2 = np.array([np.exp((self.aco2[i]) * ((self.uco2 * x) ** (self.nco2[i]))) for x in m])
        tch4 = np.array([np.exp((self.ach4[i]) * ((self.uch4 * x) ** (self.nch4[i]))) for x in m])
        tno2 = np.array([np.exp((self.ano2[i]) * ((self.uno2 * x) ** (self.nno2[i]))) for x in m])
        tco = np.array([np.exp((self.aco[i]) * ((self.uco * x) ** (self.nco[i]))) for x in m])
        # print(th2o , to3 , to2 , tco2 , tch4 , tco , tno2)
        return th2o * to3 * to2 * tco2 * tch4 * tco * tno2


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
