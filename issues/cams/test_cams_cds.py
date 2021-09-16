import os

import numpy as np
import pandas as pd
import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import curve_fit

cams_file = '/datalake/watcal/ECMWF/CAMS/2021/2021-03_month_cams-global-atmospheric-composition-forecasts.nc'
lonmin, lonmax, latmin, latmax = 113.9, 114.9, 21.6, 22.6



class lut:
    '''Load LUT FROM RT COMPUTATION (OSOAA_h)
    '''

    def __init__(self, wls=[442.3110, 492.1326, 558.9499, 664.9380, 703.8308, 739.1290,
                                    779.7236, 832.9462, 863.9796, 1610.4191, 2185.6988]):

        self.smac_bands = []
        self.N = len(wls)
        N = self.N
        self.lut_genearator = "OSOAA_h"
        self.wl = wls
        self.Cext = []
        self.Cext550 = 0
        self.vza = []
        self.sza = []
        self.azi = []
        self.aot = []
        self.refl = []

    def load_lut(self, lut_file, ind_wl=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], aot=[0.01, 0.05, 0.1, 0.3, 0.5, 0.8], vza_max=20, reflectance=True):
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
                                array of dims: [wl, sza, azi, vza]
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
            refl = nrad
            for i in range(len(self.sza)):
                refl[:, :, i, :, :] = nrad[:, :, i, :, :] / np.cos(np.radians(self.sza[i]))
            nrad = refl

        # reshape lut array for each wavelength
        N = range(len(ind_wl))
        self.refl = [[] for i in N]
        for i in N:
            self.refl[i] = nrad[:, i, :, :, :]

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

    def func(self, wl, a, b, c):
        '''function for spectral variation of AOT'''
        lnwl = np.log(wl)
        return np.exp(a + b * lnwl + c * lnwl ** 2)

    def fit_spectral_aot(self, wl, aot):
        '''call to get fitting results on AOT data'''
        self.popt, pcov = curve_fit(self.func, wl, aot)

    def get_spectral_aot(self, wl):
        '''set aot for a given set of wavelengths'''
        return self.func(wl, *self.popt)

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


def subset_xr(ds, lonmin, lonmax, latmin, latmax, lat='latitude', lon='longitude'):
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

lut_root=os.path.abspath('/work/ALT/swot/aval/OBS2CO/git/grs2/grsdata/LUT')
aero = 'rg0.10_sig0.46'
lutfine = os.path.join(lut_root,
                            'S2B/lut_' + 'osoaa_band_integrated_aot0.01_aero_' + aero + '_ws2_pressure1015.2.nc')
aero = 'rg0.80_sig0.60'
lutcoarse = os.path.join(lut_root,
                              'S2B/lut_' + 'osoaa_band_integrated_aot0.01_aero_' + aero + '_ws2_pressure1015.2.nc')

lutf = lut()
lutc = lut()
lutf.load_lut(lutfine)
lutc.load_lut(lutcoarse)
wlsat=lutf.wl
date = dt.datetime(2021, 3, 25, 2, 55)
day = date.strftime(date.strftime('%Y-%m-%d'))
wls=[380, 400,440,500,550,645,670,800,865,1020,1240,1640,2130]
idx550=4
N = len(wls)
aot = np.zeros(N, dtype=float)
ssa = np.zeros(N, dtype=float)

params=[]
for wl in wls:
    wl_=str(wl)
    params.append('aod'+wl_)
    params.append('ssa'+wl_)
cams_xr = xr.open_dataset(cams_file)[params]
cams_daily = cams_xr.sel(time=day)

#cams_sub = subset_xr(cams_daily, lonmin-1, lonmax+1, latmin-1, latmax+1)
cams_sub = subset_xr(cams_daily, lonmin, lonmax, latmin, latmax)

cams_grs = cams_sub.interp(time=date,kwargs={"fill_value": "extrapolate"})
for i, wl in enumerate(wls):
    param = 'ssa' + str(wl)
    ssa[i] = cams_grs[param].mean().data
    param = 'aod'+str(wl)
    aot[i]= cams_grs[param].mean().data

aot_sca = ssa*aot
aot550=aot[idx550]
aero = aerosol()
aero_sca = aerosol()

aero_sca.fit_spectral_aot(wls,aot_sca)
aot_grs = aero_sca.get_spectral_aot
aot_grs550=aot_grs(550)
aero.fit_spectral_aot(wls,aot)
aot_tot = aero.get_spectral_aot
aot_tot550=aot_tot(550)

# normalization of Cext to get spectral dependence of fine and coarse modes
nCext_f = lutf.Cext / lutf.Cext550
nCext_c = lutc.Cext / lutc.Cext550
print('param aerosol', nCext_f, nCext_c, aot)

aero_sca.fit_aero(nCext_f, nCext_c, aot_grs(wlsat) / aot_grs550)
aero_sca.fcoef
aero.fit_aero(nCext_f, nCext_c, aot_tot(wlsat) / aot_tot550)
aero.fcoef
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(wls,aot/aot550,'--ok')
ax1.plot(wlsat,aot_tot(wlsat) / aot_tot550,'--og')
ax1.plot(wlsat,nCext_f,':ob')
ax1.plot(wlsat,nCext_c,':or')


wl_ = np.linspace(400,2500,1001)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.plot(wls,aot,'--o')
ax1.plot(wlsat,aot_grs(wlsat),'--og')
ax1.plot(wl_,aot_tot(wl_),'-r')
ax1.plot(wl_,aot_grs(wl_),'-g')

#ax1.plot(wls,ssa,'-or')

w,h = 512,512
plt.figure()
cams_grs.aod500.plot(cmap=plt.cm.Spectral_r,vmin=0.4,vmax=1)
plt.figure()
cams_grs.aod500.interp(longitude=np.linspace(lonmin, lonmax, w),
                       latitude=np.linspace(latmax, latmin, h),
                       kwargs={"fill_value": "extrapolate"}).plot(cmap=plt.cm.Spectral_r,vmin=0.4,vmax=1)
aero = acutils.aerosol()