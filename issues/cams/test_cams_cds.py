import os

import numpy as np
import pandas as pd
import dask
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

    def load_lut(self, lut_file, ind_wl=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], aot=[0.01, 0.05, 0.1, 0.3, 0.5, 0.8],
                 vza_max=20, reflectance=True):
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

    def func(self, lnwl, a, b, c):
        '''function for spectral variation of AOT'''
        #lnwl = np.log(wl)
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


# hongkong
lonmin, lonmax, latmin, latmax = 113.9, 114.9, 21.6, 22.6
date = dt.datetime(2021, 3, 20, 10, 35)
#serponcon
lonmin, lonmax, latmin, latmax = 6.2,6.5, 44.4, 44.6
date = dt.datetime(2021, 3, 20, 10, 35)

lonmin, lonmax, latmin, latmax = -0.5,.5, 44.4, 44.6
date = dt.datetime(2021, 3, 20, 10, 35)

#22KGV

lonmin, lonmax, latmin, latmax = -49,-47, -23.57, -22.6
date = dt.datetime(2021, 3, 15, 13, 22)


year = str(date.year)
month = str(date.month).zfill(2)
cams_file = '/datalake/watcal/ECMWF/CAMS/'+year+'/'+year+'-'+month+'_month_cams-global-atmospheric-composition-forecasts.nc'

lut_root = os.path.abspath('/work/ALT/swot/aval/OBS2CO/git/grs2/grsdata/LUT')
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
wlsat = lutf.wl


day = date.strftime(date.strftime('%Y-%m-%d'))
wls = [ 400, 440, 500, 550, 645, 670, 800, 865, 1020, 1240, 1640, 2130]
rot=np.array([0.23745233, 0.15521662, 0.09104522, 0.04485828, 0.03557663,
              0.02918357, 0.02351827, 0.01829234, 0.01554304, 0.00127606, 0.00037652])
idx550 = 4
N = len(wlsat)

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

lonmin,lonmax=lonmin%360,lonmax%360
if (lonmin > 350) and (lonmax < 10):
    lonmin = lonmin - 360
    cams_daily = cams_daily.assign_coords(longitude=(((cams_daily.longitude + 180) % 360) - 180)).sortby('longitude')

lonslats = (lonmin, lonmax, latmin, latmax)

# cams_sub = subset_xr(cams_daily, lonmin-1, lonmax+1, latmin-1, latmax+1)
#cams_sub = subset_xr(cams_daily, lonmin, lonmax, latmin, latmax)
cams_sub = cams_daily.interp(longitude=np.linspace(lonmin, lonmax, 12),
                       latitude=np.linspace(latmin, latmax, 12),
                       kwargs={"fill_value": "extrapolate"})
cams_sub.aod550.plot(col='time',col_wrap=4)

# mask_lon = (cams_daily.longitude >= lonmin) & (cams_daily.longitude <= lonmax)
# mask_lat = (cams_daily.latitude >= latmin) & (cams_daily.latitude <= latmax)
# cams_ =cams_daily.where(mask_lon & mask_lat, drop=True)
# cams_.aod550.plot(col='time',col_wrap=4)

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
wl_cams = cams_aod.wavelength.str.replace('aod', '').astype(np.float64)
cams_aod = cams_aod.assign_coords(wavelength=wl_cams)

aero = aerosol()

r, c = cams_grs.aod550.shape
aot_550 = np.zeros((r, c), dtype=float)
fcoef = np.zeros((r, c), dtype=float)
aot_ = np.zeros((N,r, c), dtype=float)
ssa_grs = np.zeros((r, c, N), dtype=float)
for ir in range(r):
    for ic in range(c):
        print(ir)
        aero.fit_spectral_aot(cams_aod.wavelength, cams_aod.data[...,ir,ic])
        aot_[:,ir,ic] = aero.get_spectral_aot(wlsat)

ssa_grs = cams_ssa.interp(wavelength=wlsat,kwargs={"fill_value": "extrapolate"}).chunk({'wavelength':1})
aot_grs = ssa_grs.copy(data=aot_).chunk({'wavelength':1})
aot_sca_grs = (ssa_grs * aot_grs).chunk({'wavelength':1})
aot_sca_550 = aot_sca_grs.interp(wavelength=550,method='cubic')

ssa = np.array((aot_sca_grs[:,0,0]+rot)/(aot_grs[:,0,0]+rot))

# normalization of Cext to get spectral dependence of fine and coarse modes
nCext_f = lutf.Cext / lutf.Cext550
nCext_c = lutc.Cext / lutc.Cext550
for ir in range(r):
    for ic in range(c):
        fcoef[ir,ic] = aero.fit_aero(nCext_f, nCext_c, aot_sca_grs[...,ir,ic] / aot_sca_550[ir,ic])[0]
fcoef = aot_sca_550.copy(data=fcoef)


w, h = 549, 549
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.subplots_adjust(bottom=0.1, top=0.95, left=0.1, right=0.99,
                    hspace=0.15, wspace=0.1)
fcoef.plot(ax=axs[0,0],cmap=plt.cm.Spectral_r)
fcoef.interp(longitude=np.linspace(lonmin, lonmax, 12),
                       latitude=np.linspace(latmax, latmin, 12),
                       kwargs={"fill_value": "extrapolate"}).interp(longitude=np.linspace(lonmin, lonmax, 512),
                       latitude=np.linspace(latmax, latmin, 512),
                       kwargs={"fill_value": "extrapolate"}).interp(longitude=np.linspace(lonmin, lonmax, w),
                       latitude=np.linspace(latmax, latmin, h),method="nearest",
                       kwargs={"fill_value": "extrapolate"}).plot(ax=axs[0,1],cmap=plt.cm.Spectral_r)
aot_sca_grs.sel(wavelength=560).plot(ax=axs[1,0],cmap=plt.cm.Spectral_r)
axs[1,0].set_title('AOT_sca')
aot_sca_grs.sel(wavelength=560).interp(longitude=np.linspace(lonmin, lonmax, w),
                       latitude=np.linspace(latmax, latmin, h),
                       kwargs={"fill_value": "extrapolate"}).plot(ax=axs[1,1],cmap=plt.cm.Spectral_r)
axs[1,1].set_title('interp. AOT_sca')

aot_grs.sel(wavelength=560).interp(longitude=np.linspace(lonmin, lonmax, w),
                       latitude=np.linspace(latmax, latmin, h),
                       kwargs={"fill_value": "extrapolate"}).plot(ax=axs[2,0],cmap=plt.cm.Spectral_r, vmin=0.4, vmax=0.61)
axs[2,0].set_title('interp. AOT_sca')

ssa_grs.sel(wavelength=560).interp(longitude=np.linspace(lonmin, lonmax, w),
                       latitude=np.linspace(latmax, latmin, h),
                       kwargs={"fill_value": "extrapolate"}).plot(ax=axs[2,1],cmap=plt.cm.Spectral_r) #, vmin=0.8, vmax=1)
axs[2,1].set_title('interp. ssa')

################################
#########""
###############################
ir,ic=0,0
aot=cams_aod.data[:,ir,ic]
aero.fit_spectral_aot(wls, aot)
aot_grs = aero.get_spectral_aot
aot_grs550 = aot_grs(550)


# normalization of Cext to get spectral dependence of fine and coarse modes
nCext_f = lutf.Cext / lutf.Cext550
nCext_c = lutc.Cext / lutc.Cext550
print('param aerosol', nCext_f, nCext_c, aot)

aero.fit_aero(nCext_f, nCext_c, aot_grs(wlsat) / aot_grs550)
aero.fcoef


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
axs[0].plot(wlsat, aot_grs(wlsat) / aot_grs550, '--ok',label='CAMS')
axs[0].plot(wlsat, nCext_f, ':ob',label='fine')
axs[0].plot(wlsat, nCext_c, ':or',label='coarse')
axs[0].plot(wlsat, aero.fcoef * nCext_f + (1 - aero.fcoef) * nCext_c, '*',label='fcoef:'+'{:5.3f}'.format(aero.fcoef[0]))
axs[0].set_xlabel('Wavelength (nm)')
axs[0].set_ylabel('nCext or naot')

axs[0].legend()

axs[1].plot(wlsat, lutf.Cext, ':ob',label='Cext fine')
axs[1].plot(wlsat, lutc.Cext, ':or',label='Cext coarse')
axs[1].plot(wlsat, aero.fcoef * lutf.Cext + (1 - aero.fcoef) * lutc.Cext, ':*',label='Cext mixture')
a_=axs[1].twinx()
a_.plot(wlsat, aot_grs(wlsat), '--ok',label='aot')
axs[1].set_xlabel('Wavelength (nm)')
axs[1].set_ylabel('Cext or ssa')
axs[1].legend()

cams_aod.plot(col='wavelength',col_wrap=4)
cams_ssa.plot(col='wavelength',col_wrap=4)

Naot=lutc.aot.__len__()
rtoaf = np.zeros((Naot, N))
rtoac = np.zeros((Naot, N))
for iband in range(N):
    for iaot in range(Naot):
        rtoaf[iaot, iband] = lutf.refl[iband][iaot, 0,0,1]
        rtoac[iaot, iband] = lutc.refl[iband][iaot, 0,0,1]
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
iaot=2
rtoa_tot = aero.fcoef *rtoaf[iaot]+(1-aero.fcoef) *rtoac[iaot]
for i in [0,1]:
    axs[i].plot(wlsat,rtoaf[iaot],'o--b',label='fine')
    axs[i].plot(wlsat,rtoac[iaot],'o--r',label='coarse')
    axs[i].plot(wlsat,rtoa_tot,'o-k',label='tot')
    axs[i].plot(wlsat,ssa_grs[:,ir,ic]*rtoa_tot,'o:g',label='ssa adjusted')
    axs[i].legend()
axs[1].semilogy()
