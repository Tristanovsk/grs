import os

import importlib_resources
import yaml

import numpy as np
import xarray as xr

# keep attributes through operation on xarray objects
xr.set_options(keep_attrs=True)
import rioxarray as rio
import logging
import gc

from numba import njit
from scipy import ndimage

from multiprocessing import Pool  # Process pool
from multiprocessing import sharedctypes
import itertools

from . import AuxData, Rasterization, GaseousTransmittance, Misc


class Kernel():
    def __init__(self,
                 prod,
                 aero_lut,
                 trans_lut,
                 cams):

        self.prod = prod
        self.raster = prod.raster
        self.aero_lut = aero_lut
        self.trans_lut = trans_lut
        self.aot_ref_lim = np.max(aero_lut.aot_ref.values)

        self.cams = cams
        aot_ref_cams = cams.cams_aod.interp(wl=550)
        self.aot_ref_cams_max = aot_ref_cams.max()
        self.aot_ref_cams_mean = aot_ref_cams.mean()
        self.aot_ref_cams_min = aot_ref_cams.min()

        self.wl_true = prod.raster.wl_true
        self.pressure_ref = 101320.

        self.neg_pix_max = 0

        if 'S2' in prod.sensor:
            self.monoview = False
        else:
            self.monoview = True
        self._R_ = Rasterization(monoview=self.monoview)

    def get_coarse_masked_raster(self,
                                 variables=['sza', 'vza', 'raa', 'bands'],

                                 mask=None,
                                 xcoarsen=2,
                                 ycoarsen=2,
                                 min_pixel_prop=0.5):
        if mask is not None:
            raster = self.raster[variables].where(mask == 0)
        else:
            raster = self.raster[variables]
        self.coarse_masked_raster = raster.coarsen(x=xcoarsen, y=ycoarsen, boundary="pad").mean()
        self.coarse_masked_raster['water_pixel_number'] = raster['sza']. \
            coarsen(x=xcoarsen, y=ycoarsen, boundary="pad").count()
        min_pixel_number = min_pixel_prop * np.max(self.coarse_masked_raster['water_pixel_number'])
        return self.coarse_masked_raster.where(self.coarse_masked_raster.water_pixel_number > min_pixel_number)

    def lut_preparation(self,
                        wind=2,
                        ang_resol={'sza': 1, 'vza': 1, 'raa': 0},
                        aot_refs=np.linspace(0.0, 0.8, 25),
                        weights=[0, 0.5, 1, 0., 0.]):

        logging.info('LUT preparation')

        aero_lut = self.aero_lut.sel(wind=wind, method='nearest')
        trans_lut = self.trans_lut.sel(wind=wind, method='nearest')
        for param in ['sza', 'vza', 'raa']:
            self.prod.raster[param + '_trunc'] = self.prod.raster[param].round(ang_resol[param])

        # set mixture of aerosol models
        # ['ARCT_rh70', 'COAV_rh70', 'DESE_rh70', 'MACL_rh70', 'URBA_rh70']
        for ii, weight in enumerate(weights):
            if ii == 0:
                mix_aero_lut_ = weight * aero_lut.isel(model=ii)
                mix_trans_lut_ = weight * trans_lut.isel(model=ii)
            else:
                mix_aero_lut_ = mix_aero_lut_ + weight * aero_lut.isel(model=ii)
                mix_trans_lut_ = mix_trans_lut_ + weight * trans_lut.isel(model=ii)
        mix_aero_lut_ = mix_aero_lut_ / np.sum(weights)
        mix_trans_lut_ = mix_trans_lut_ / np.sum(weights)

        aero_lut = mix_aero_lut_.interp(wl=self.wl_true, method='quadratic')
        self.trans_aero_lut = mix_trans_lut_.interp(wl=self.wl_true, method='quadratic')
        del mix_aero_lut_, mix_trans_lut_

        sza_ = np.unique(self.prod.raster.sza_trunc)
        vza_ = np.unique(self.prod.raster.vza_trunc)
        azi_ = np.unique((180. - self.prod.raster.raa_trunc) % 360)

        sza_ = sza_[~np.isnan(sza_)]
        vza_ = vza_[~np.isnan(vza_)]
        azi_ = azi_[~np.isnan(azi_)]

        self.trans_aero_lut = self.trans_aero_lut.interp(sza=[*sza_, *vza_]).interp(aot_ref=aot_refs,
                                                                                    method='quadratic')

        Rdiff_lut = aero_lut.I.interp(sza=sza_, vza=vza_).interp(azi=azi_).interp(aot_ref=aot_refs, method='quadratic')
        self.aot_lut = aero_lut.aot.interp(aot_ref=aot_refs, method='quadratic')
        self.Rdiff_lut = Rdiff_lut
        self.Rray = Rdiff_lut.sel(aot_ref=0)

        self.szas = Rdiff_lut.sza.values
        self.vzas = Rdiff_lut.vza.values
        self.azis = Rdiff_lut.azi.values
        self.aot_refs = Rdiff_lut.aot_ref.values

        _auxdata = AuxData(wl=self.wl_true)  # wl=masked.wl)
        self.sunglint_eps = _auxdata.sunglint_eps  # ['mean'].interp(wl=wl_true)
        self.rot = _auxdata.rot

    def mu_N(self,
             sza,
             vza,
             azi,
             monoview=False):
        '''
        Compute the normal angle to wave slopes that produce sunglint.
        Warning: azi: azimuth in rad for convention azi=180 when sun-sensenor in oppositioon

        :param sza: solar zenith angle in degree
        :param vza: viewing zenith angle in degree
        :param azi: relative azimuth bewteen sun and sensor
        :param monoview: Set sensor viewing  configuration:
            - monoview = True : same viewing angles for all the spectral bands
            - monoview = False : viewing angles depend on spectral band (e.g. Sentinel-2 images)

        :return: cosine of normal angle to sunglint wave facets
        '''
        vzar = np.radians(vza)
        azir = np.radians(azi)
        szar = np.radians(sza)
        cos_alpha = np.cos(azir) * np.sin(vzar) * np.sin(szar) + np.cos(vzar) * np.cos(szar)
        xmu_N = (np.cos(szar) + np.cos(vzar)) / np.sqrt(2 * (1 + cos_alpha))
        if monoview:
            xmu_N = xmu_N.transpose('y', 'x')
        else:
            # ensure similar shape as inputs
            xmu_N = xmu_N.transpose('wl', 'y', 'x')
        return xmu_N

    def p_slope(self,
                sza,
                vza,
                azi,
                sigma2=0.02,
                monoview=False):
        '''
        Compute propability of wave slopes producing sunglint.

        :param sza: solar zenith angle in degree
        :param vza: viewing zenith angle in degree
        :param azi: relative azimuth bewteen sun and sensor
        :param sigma2: mean square slope of the wave slope distribution
        :param monoview: Set sensor viewing  configuration:
            - monoview = True : same viewing angles for all the spectral bands
            - monoview = False : viewing angles depend on spectral band (e.g. Sentinel-2 images)

        :return:
        '''

        cosN = self.mu_N(sza, vza, azi, monoview=monoview)
        thetaN = np.arccos(cosN)
        # stats == 'cm_iso':
        # TODO check consitency between sigma2 and formulation
        # Pdist_ = 1. / (np.pi *2.* sigma2) * np.exp(-1./2 * np.tan(thetaN) ** 2 / sigma2)
        xp_slope = 1. / (np.pi * sigma2) * np.exp(- np.tan(thetaN) ** 2 / sigma2) / cosN ** 4
        if monoview:
            xp_slope = xp_slope.transpose('y', 'x')
        else:
            xp_slope = xp_slope.transpose('wl', 'y', 'x')
        return xp_slope

    def set_gas_transmittance(self):
        logging.info('set gaseous transmittance')
        # -------------------------------------------------------------
        # SET GASEOUS TRANSMITTANCE FOR LOW ALTITUDE GASES
        # -------------------------------------------------------------
        gases = ['h2o', 'ch4']
        gas_trans = GaseousTransmittance(self.prod, self.cams)
        # set total transmittance values
        gas_trans.coef_abs_scat['h2o'] = 0.5
        gas_trans.coef_abs_scat['ch4'] = 0.65
        self.Tg_diff_raster = gas_trans.get_gaseous_transmittance(gases=gases, background=False).transpose("wl", "y",
                                                                                                           "x")

        # set total transmittance values
        for gas in gases:
            gas_trans.coef_abs_scat[gas] = 1
        self.Tg_raster = gas_trans.get_gaseous_transmittance(gases=gases, background=False).transpose("wl", "y", "x")

    def rugosity_est_chunk(self,
                           raster,
                           aot_ref=0.1,
                           winds=np.linspace(0.2, 15, 25),
                           idx_nir=[6, 7, 8, 10],
                           chunk=100000
                           ):

        monoview = self.monoview
        _R_ = self._R_

        raster = raster.isel(wl=idx_nir)
        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        if Ny == 0 or Nx == 0:
            return

        _wl_true = self.wl_true[idx_nir]
        Rdiff_lut_ = self.Rdiff_lut.isel(wl=idx_nir)

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.isel(wl=idx_nir).interp(x=raster.x, y=raster.y)
        _Tg_raster = self.Tg_raster.isel(wl=idx_nir).interp(x=raster.x, y=raster.y)

        # subsetting
        _sza = raster.sza
        _raa = raster.raa
        _vza = raster.vza
        _vza_mean = np.mean(_vza, axis=0).values
        _azi = (180. - _raa) % 360

        _pressure_ = _pressure / self.pressure_ref

        _Tg_abs = _Tg_raster.values
        _Tg_abs_diff = _Tg_diff_raster.values

        # get LUT values
        _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza.values,
                                          self.vzas, _vza.values,
                                          self.azis, _azi.values,
                                          Nwl, Ny, Nx,
                                          Rdiff_lut_.interp(aot_ref=aot_ref, method='linear').values)

        _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

        #  correction for diffuse light
        Rcorr = _band_rad - _Rdiff

        delta_R = []
        for wind in winds:
            _sigma2 = (wind + 0.586) / 195.3

            _p_slope = self.p_slope(_sza, _vza, _raa, sigma2=_sigma2,
                                    monoview=monoview).values  # _p_slope[:, iy:yc,ix:xc]

            Rgeom = _p_slope[2] * Rcorr / (_Tg_abs * _p_slope)

            deltaR = ((Rgeom[2] - Rgeom[0]) * (_wl_true[1] - _wl_true[0]) /
                      (_wl_true[2] - _wl_true[0]) + Rgeom[0] - Rgeom[1])
            delta_R.append(deltaR.assign_coords({'wind': wind}))
        swir_Rcorr = Rcorr[-1]

        delta_R = xr.concat(delta_R, dim='wind')

        @njit
        def find_root(cost, nx, ny, values, result):
            for ix in range(nx):
                for iy in range(ny):

                    if np.isnan(cost[:, iy, ix]).all():
                        continue

                    ind_wind = np.argmin(cost[:, iy, ix])

                    result[ix, iy, :] = [values[ind_wind], cost[ind_wind, iy, ix]]
            return result

        cost = delta_R.values ** 2

        values = winds
        result = np.full((Nx, Ny, 2), np.nan)
        result = find_root(cost, Nx, Ny, values, result)

        return xr.Dataset(dict(wind=(["y", "x"], result[:, :, 0].T),
                               cost=(["y", "x"], result[:, :, 1].T),
                               brdf_g=(["y", "x"], swir_Rcorr.values)
                               ),
                          coords=dict(
                              x=raster.x,
                              y=raster.y),
                          attrs=dict(
                              description="wind speed retrieved from sunglint, parallax and isotropic Cox-Munk model",
                              units='m.s-1')
                          )

    def rugosity_est(self,
                     raster,
                     aot_ref=0.1,
                     winds=np.linspace(0.2, 15, 25),
                     idx_nir=[6, 7, 8, 10],
                     chunk=100000
                     ):

        monoview = self.monoview
        _R_ = self._R_

        _raster = raster.isel(wl=idx_nir)
        _Nwl, _height, _width = _raster.bands.shape
        _wl_true = self.wl_true[idx_nir]
        Rdiff_lut_ = self.Rdiff_lut.isel(wl=idx_nir)

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.interp(x=raster.x, y=raster.y)
        _Tg_raster = self.Tg_raster.interp(x=raster.x, y=raster.y)

        # TODO check if possible/necessary to implement the chunk process (not working now)
        for iy in range(0, _height, chunk):
            yc = min(_height, iy + chunk)

            for ix in range(0, _width, chunk):

                xc = min(_width, ix + chunk)

                _band_rad = _raster.bands[:, iy:yc, ix:xc]

                Nwl, Ny, Nx = _band_rad.shape
                if Ny == 0 or Nx == 0:
                    continue
                arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

                # subsetting
                _sza = _raster.sza[iy:yc, ix:xc]  # .values
                _raa = _raster.raa[:, iy:yc, ix:xc]
                _vza = _raster.vza[:, iy:yc, ix:xc]
                _vza_mean = np.mean(_vza, axis=0).values
                _azi = (180. - _raa) % 360

                _pressure_ = _pressure[iy:yc, ix:xc] / self.pressure_ref

                _Tg_abs = _Tg_raster[idx_nir, iy:yc, ix:xc].values
                _Tg_abs_diff = _Tg_diff_raster[idx_nir, iy:yc, ix:xc].values

                # get LUT values
                _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza.values,
                                                  self.vzas, _vza.values,
                                                  self.azis, _azi.values,
                                                  Nwl, Ny, Nx,
                                                  Rdiff_lut_.interp(aot_ref=aot_ref, method='linear').values)

                _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

                #  correction for diffuse light
                Rcorr = _band_rad - _Rdiff

                delta_R = []
                for wind in winds:
                    _sigma2 = (wind + 0.586) / 195.3

                    _p_slope = self.p_slope(_sza, _vza, _raa, sigma2=_sigma2,
                                            monoview=monoview).values  # _p_slope[:, iy:yc,ix:xc]

                    Rgeom = _p_slope[2] * Rcorr / (_Tg_abs * _p_slope)

                    deltaR = ((Rgeom[2] - Rgeom[0]) * (_wl_true[1] - _wl_true[0]) /
                              (_wl_true[2] - _wl_true[0]) + Rgeom[0] - Rgeom[1])
                    delta_R.append(deltaR.assign_coords({'wind': wind}))
        swir_Rcorr = Rcorr[-1]

        delta_R = xr.concat(delta_R, dim='wind')

        @njit
        def find_root(cost, nx, ny, values, result):
            for ix in range(nx):
                for iy in range(ny):

                    if np.isnan(cost[:, iy, ix]).all():
                        continue

                    ind_wind = np.argmin(cost[:, iy, ix])

                    result[ix, iy, :] = [values[ind_wind], cost[ind_wind, iy, ix]]
            return result

        cost = delta_R.values ** 2
        nx, ny = _width, _height
        values = winds
        result = np.full((_width, _height, 2), np.nan)
        result = find_root(cost, nx, ny, values, result)

        self.wind_img = xr.Dataset(dict(wind=(["y", "x"], result[:, :, 0].T),
                                        cost=(["y", "x"], result[:, :, 1].T),
                                        brdf_g=(["y", "x"], swir_Rcorr.values)
                                        ),
                                   coords=dict(
                                       x=raster.x,
                                       y=raster.y),
                                   attrs=dict(
                                       description="wind speed retrieved from sunglint, parallax and isotropic Cox-Munk model",
                                       units='m.s-1')
                                   )

    def aerosol_swir_chunk(self,
                           raster,
                           aot_refs=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8],
                           idx_swir=[-2, -1],
                           neg_pix_rel_thresh=0.05,
                           chunk=100000
                           ):

        monoview = self.monoview
        _R_ = self._R_

        raster = raster.isel(wl=idx_swir)
        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        if Ny == 0 or Nx == 0:
            return

        _wl_true = self.wl_true[idx_swir]
        Rdiff_lut_ = self.Rdiff_lut.isel(wl=idx_swir)

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.isel(wl=idx_swir).interp(x=raster.x, y=raster.y)
        _Tg_abs_diff = _Tg_diff_raster.values

        sigma2 = (self.wind_img.wind + 0.586) / 195.3
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)

        # subsetting
        _sza = raster.sza
        _raa = raster.raa
        _vza = raster.vza
        _vza_mean = np.mean(_vza, axis=0).values
        _azi = (180. - _raa) % 360

        _pressure_ = _pressure / self.pressure_ref

        xres = []
        for _aot_ref in aot_refs:
            _aot = self.aot_lut.interp(aot_ref=_aot_ref, method='quadratic')  # .values

            # get LUT values
            _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza.values,
                                              self.vzas, _vza.values,
                                              self.azis, _azi.values,
                                              Nwl, Ny, Nx, Rdiff_lut_.sel(aot_ref=_aot_ref, method='nearest').values)

            _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

            #  correction for diffuse light
            Rcorr = _band_rad.values - _Rdiff

            xres_ = xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rcorr),
                                    ),
                               coords=dict(wl=raster.wl,
                                           x=raster.x,
                                           y=raster.y,
                                           aot_ref=_aot_ref), ).__deepcopy__()
            xres.append(xres_)

        self.xres = xr.concat(xres, dim='aot_ref')

        # get aot550 max from which negative Rrs appear
        metric = self.xres.Rrs.min('wl')
        metric = metric.where(metric.isel(aot_ref=0) > 0)

        metric_pixnum = metric.where(metric < 0).count(['x', 'y'])
        neg_pix_max = metric_pixnum.max()
        metric_pixnum_max = neg_pix_max * neg_pix_rel_thresh

        aot_ref_maxs = metric_pixnum.where(metric_pixnum < metric_pixnum_max, drop=True).aot_ref
        if len(aot_ref_maxs) > 0:
            self.aot_ref_max = aot_ref_maxs.values[-1]
        else:
            self.aot_ref_max = aot_refs[-1]*0.8 #self.aot_ref_cams_max

        # for further lut interpolation:
        #self.aot_ref_max = np.max([0.01, self.aot_ref_max])

        if ~((self.aot_ref_max > 0) & (self.aot_ref_max < self.aot_ref_lim)):
            self.aot_ref_max = np.min([self.aot_ref_cams_max, self.aot_ref_lim])

    def aerosol_swir(self,
                     raster,
                     aot_refs=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8],
                     idx_swir=[-2, -1],
                     neg_pix_rel_thresh=0.05,
                     chunk=100000
                     ):

        monoview = self.monoview
        _R_ = self._R_

        _raster = raster.isel(wl=idx_swir)
        _Nwl, _height, _width = _raster.bands.shape
        Rdiff_lut_ = self.Rdiff_lut.isel(wl=idx_swir).interp(aot_ref=aot_refs, method='quadratic')

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.interp(x=raster.x, y=raster.y)

        sigma2 = (self.wind_img.wind + 0.586) / 195.3
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)

        Rrs_tmp = np.full((_Nwl, _height, _width), np.nan, dtype=self.prod._type)

        xres = []
        for _aot_ref in aot_refs:
            _aot = self.aot_lut.interp(aot_ref=_aot_ref, method='quadratic')  # .values
            for iy in range(0, _height, chunk):
                yc = min(_height, iy + chunk)

                for ix in range(0, _width, chunk):
                    xc = min(_width, ix + chunk)

                    _band_rad = _raster.bands[:, iy:yc, ix:xc]

                    Nwl, Ny, Nx = _band_rad.shape
                    if Ny == 0 or Nx == 0:
                        continue
                    arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

                    # subsetting
                    _sigma2_ = _sigma2[iy:yc, ix:xc]
                    _sza = _raster.sza[iy:yc, ix:xc]  # .values
                    if monoview:
                        _raa = _raster.raa[iy:yc, ix:xc]
                        _vza = _raster.vza[iy:yc, ix:xc]
                        _vza_mean = _vza.values
                    else:
                        _raa = _raster.raa[:, iy:yc, ix:xc]
                        _vza = _raster.vza[:, iy:yc, ix:xc]
                        _vza_mean = np.mean(_vza, axis=0).values

                    _azi = (180. - _raa) % 360

                    _pressure_ = _pressure[iy:yc, ix:xc] / self.pressure_ref

                    _Tg_abs_diff = _Tg_diff_raster[idx_swir, iy:yc, ix:xc].values

                    # get LUT values
                    _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza.values,
                                                      self.vzas, _vza.values,
                                                      self.azis, _azi.values,
                                                      Nwl, Ny, Nx,
                                                      Rdiff_lut_.sel(aot_ref=_aot_ref, method='nearest').values)

                    _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

                    #  correction for diffuse light
                    Rcorr = _band_rad.values - _Rdiff
                    Rrs_tmp[:, iy:yc, ix:xc] = Rcorr

            xres_ = xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rrs_tmp),
                                    ),
                               coords=dict(wl=_raster.wl,
                                           x=_raster.x,
                                           y=_raster.y,
                                           aot_ref=_aot_ref), ).__deepcopy__()
            xres.append(xres_)

        self.xres = xr.concat(xres, dim='aot_ref')

        # get aot550 max from which negative Rrs appear
        metric = self.xres.Rrs.min('wl')
        metric = metric.where(metric.isel(aot_ref=0) > 0)

        metric_pixnum = metric.where(metric < 0).count(['x', 'y'])
        neg_pix_max = metric_pixnum.max()
        metric_pixnum_max = neg_pix_max * neg_pix_rel_thresh

        aot_ref_maxs = metric_pixnum.where(metric_pixnum < metric_pixnum_max, drop=True).aot_ref
        if len(aot_ref_maxs) > 0:
            self.aot_ref_max = aot_ref_maxs.values[-1]
        else:
            self.aot_ref_max = aot_refs[-1] #self.aot_ref_cams_max

        # for further lut interpolation:
        # self.aot_ref_max = np.max([0.011, aot_ref_max])

        # to stay within LUT aot range
        if ~((self.aot_ref_max > 0) & (self.aot_ref_max < self.aot_ref_lim)):
            self.aot_ref_max = np.min([self.aot_ref_cams_max, self.aot_ref_lim])

        # TODO uncomment after end of dev
        # del self.xres

    def aerosol_visible_chunk(self,
                              raster,
                              # aot_refs=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8],
                              iwl_swir=[-2, -1],
                              neg_pix_rel_thresh=0.05,
                              chunk=100000
                              ):
        ##############################################
        # Estimate AOT from visible-NIR
        ##############################################
        monoview = self.monoview
        _R_ = self._R_

        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        if Ny == 0 or Nx == 0:
            return

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.interp(x=raster.x, y=raster.y)
        _Tg_abs_diff = _Tg_diff_raster.values
        _Tg_raster = self.Tg_raster.interp(x=raster.x, y=raster.y)
        _Tg_abs = _Tg_raster.values

        # parameters
        _sza = raster.sza
        _raa = raster.raa
        _vza = raster.vza
        _vza_mean = np.mean(_vza, axis=0).values
        _azi = (180. - _raa) % 360
        _air_mass_ = Misc.air_mass(_sza, _vza)
        _pressure_ = _pressure / self.pressure_ref

        sigma2 = (self.wind_img.wind + 0.586) / 195.3
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)
        _p_slope_ = self.p_slope(_sza, _vza, _raa, sigma2=_sigma2,
                                 monoview=monoview).values

        xres = []

        sunglint_eps = self.sunglint_eps.values

        aot_ref_bound = np.min([self.aot_ref_max * 1.2, self.aot_ref_lim])
        aot_refs = np.linspace(0, aot_ref_bound, 21)

        xres = []
        for _aot_ref in aot_refs:
            _aot = self.aot_lut.interp(aot_ref=_aot_ref, method='quadratic')  # .values

            arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

            # get LUT values
            _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza.values,
                                              self.vzas, _vza.values,
                                              self.azis, _azi.values,
                                              Nwl, Ny, Nx,
                                              self.Rdiff_lut.sel(aot_ref=_aot_ref, method='nearest').values)

            _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

            #  correction for diffuse light
            Rcorr = _band_rad.values - _Rdiff

            # direct transmittance up/down
            Tdir = np.exp(-(_aot + self.rot.values * np.mean(
                _pressure_)) * _air_mass_)  # acutils.Misc.transmittance_dir(_aot, _air_mass_, _rot_raster)
            self.Tdir = Tdir
            self._p_slope_ = _p_slope_
            self._Tg_abs = _Tg_abs
            self._p_slope_ = _p_slope_
            Rf = np.full((len(iwl_swir), Ny, Nx), np.nan, dtype=np.float32)
            for iwl in iwl_swir:
                if monoview:
                    Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_abs[iwl] * sunglint_eps[iwl] * _p_slope_)
                else:
                    Rf[iwl] = (sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl])
                    Rf[iwl] = Rf[iwl] / (
                            _Tg_abs[iwl] * Tdir[iwl] * sunglint_eps[iwl] * _p_slope_[iwl])

            Rf[Rf < -0.0004] = np.nan
            Rf[Rf < 0.] = 0.0
            Rf = np.min(Rf, axis=0)

            Rf = _R_._multiplicate(sunglint_eps, Rf, arr_tmp)
            Rf = _Tg_abs * Tdir * Rf * _p_slope_ / (sunglint_eps[-1] * _p_slope_[-1])

            # sunglint removal
            Rrs_toa = ((Rcorr - Rf) / np.pi)

            # print('success')
            xres_ = xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rrs_toa.values),
                                    ),
                               coords=dict(wl=raster.wl,
                                           x=raster.x,
                                           y=raster.y,
                                           aot_ref=_aot_ref),
                               ).__deepcopy__()
            xres.append(xres_)
        xres = xr.concat(xres, dim='aot_ref')
        self.xres = xres

        metric = (xres.Rrs.isel(wl=slice(0, 7)).min('wl'))
        # remove pixel with negative Rrs in Rayleigh correction (e.g., shadows)
        metric = metric.where(metric.isel(aot_ref=0) > 0)

        metric_pixnum = metric.where(metric < 0).count(['x', 'y'])
        neg_pix_max = np.max([self.neg_pix_max, metric_pixnum.max()])
        metric_pixnum_max = neg_pix_max * neg_pix_rel_thresh

        aot_ref_maxs = metric_pixnum.where(metric_pixnum < metric_pixnum_max, drop=True).aot_ref
        if len(aot_ref_maxs) > 0:
            aot_ref_max = aot_ref_maxs.values[-1]
        else:
            retrieved_aot = (self.cams.cams_aod.interp(wl=550, method='quadratic')).interp(x=raster.x, y=raster.y)
            retrieved_aot.name = 'aot_ref'
            return retrieved_aot

        ind = np.nanargmin((aot_refs - aot_ref_max) ** 2)

        map_pixels = metric.sel(aot_ref=aot_refs[ind])
        map_pixels = self.cams.cams_aod.interp(wl=550, method='quadratic').interp(x=metric.x, y=metric.y).where(
            map_pixels < 0)

        cams_aot_ref_mean_ = map_pixels.mean().values
        if np.isnan(cams_aot_ref_mean_):
            cams_aot_ref_mean_ = self.cams.cams_aod.interp(wl=550, method='quadratic').mean()

        self.recalib_aot_ref = float(aot_ref_max / cams_aot_ref_mean_)
        retrieved_aot = (self.cams.cams_aod.interp(wl=550, method='quadratic') * self.recalib_aot_ref).interp(
            x=raster.x, y=raster.y)
        retrieved_aot.name = 'aot_ref'
        return retrieved_aot

    def aerosol_visible(self,
                        raster,
                        # aot_refs=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8],
                        iwl_swir=[-2, -1],
                        neg_pix_rel_thresh=0.05,
                        chunk=100000
                        ):
        ##############################################
        # Estimate AOT from visible-NIR
        ##############################################
        monoview = self.monoview
        _R_ = self._R_

        _Nwl, _height, _width = raster.bands.shape
        sunglint_eps = self.sunglint_eps.values

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.interp(x=raster.x, y=raster.y)

        sigma2 = (self.wind_img.wind + 0.586) / 195.3
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)
        _Nwl, _height, _width = raster.bands.shape

        _pressure = self.cams.raster.sp.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.interp(x=raster.x, y=raster.y)
        _Tg_raster = self.Tg_raster.interp(x=raster.x, y=raster.y)
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)

        Rrs_tmp = np.full((_Nwl, _height, _width), np.nan, dtype=self.prod._type)
        p_slope_tmp = np.full((_Nwl, _height, _width), np.nan, dtype=self.prod._type)
        Rf_tmp = np.full((_height, _width), np.nan, dtype=self.prod._type)

        aot_ref_bound = np.min([self.aot_ref_max * 1.2, self.aot_ref_lim])
        aot_refs = np.linspace(0, aot_ref_bound, 21)

        xres = []
        for _aot_ref in aot_refs:
            _aot = self.aot_lut.interp(aot_ref=_aot_ref, method='quadratic')  # .values
            for iy in range(0, _height, chunk):
                yc = min(_height, iy + chunk)

                for ix in range(0, _width, chunk):
                    xc = min(_width, ix + chunk)

                    _band_rad = raster.bands[:, iy:yc, ix:xc]

                    Nwl, Ny, Nx = _band_rad.shape
                    if Ny == 0 or Nx == 0:
                        continue
                    arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

                    # subsetting
                    _sigma2_ = _sigma2[iy:yc, ix:xc]
                    _sza = raster.sza[iy:yc, ix:xc]  # .values
                    if monoview:
                        _raa = raster.raa[iy:yc, ix:xc]
                        _vza = raster.vza[iy:yc, ix:xc]
                        _vza_mean = _vza.values
                    else:
                        _raa = raster.raa[:, iy:yc, ix:xc]
                        _vza = raster.vza[:, iy:yc, ix:xc]
                        _vza_mean = np.mean(_vza, axis=0).values

                    _azi = (180. - _raa) % 360
                    _air_mass_ = Misc.air_mass(_sza,
                                               _vza)  # air_mass[:, iy:yc,ix:xc] #air_mass(_sza,_vza).values #_p_slope = prod.raster.p_slope[:, iy:yc,ix:xc]
                    _p_slope_ = self.p_slope(_sza, _vza, _raa, sigma2=_sigma2_,
                                             monoview=monoview).values  # _p_slope[:, iy:yc,ix:xc]

                    _pressure_ = _pressure[iy:yc, ix:xc] / self.pressure_ref
                    _Tg_abs = _Tg_raster[:, iy:yc, ix:xc].values
                    _Tg_abs_diff = _Tg_diff_raster[:, iy:yc, ix:xc].values

                    # get LUT values
                    _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza.values,
                                                      self.vzas, _vza.values,
                                                      self.azis, _azi.values,
                                                      Nwl, Ny, Nx,
                                                      self.Rdiff_lut.sel(aot_ref=_aot_ref, method='nearest').values)

                    _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

                    #  correction for diffuse light
                    Rcorr = _band_rad.values - _Rdiff

                    # direct transmittance up/down
                    Tdir = np.exp(-(_aot + self.rot.values * np.mean(
                        _pressure_)) * _air_mass_)  # acutils.Misc.transmittance_dir(_aot, _air_mass_, _rot_raster)
                    self.Tdir = Tdir
                    self._p_slope_ = _p_slope_
                    self._Tg_abs = _Tg_abs
                    self._p_slope_ = _p_slope_
                    Rf = np.full((len(iwl_swir), Ny, Nx), np.nan, dtype=np.float32)
                    for iwl in iwl_swir:
                        if monoview:
                            Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_abs[iwl] * sunglint_eps[iwl] * _p_slope_)
                        else:
                            Rf[iwl] = (sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl])
                            Rf[iwl] = Rf[iwl] / (
                                    _Tg_abs[iwl] * Tdir[iwl] * sunglint_eps[iwl] * _p_slope_[iwl])

                    Rf[Rf < -0.0004] = np.nan
                    Rf[Rf < 0.] = 0.0
                    Rf = np.min(Rf, axis=0)
                    Rf_tmp[iy:yc, ix:xc] = Rf

                    Rf = _R_._multiplicate(sunglint_eps, Rf, arr_tmp)
                    Rf = _Tg_abs * Tdir * Rf * _p_slope_ / (sunglint_eps[-1] * _p_slope_[-1])

                    # sunglint removal
                    # Rrs_tmp_ = Rrs_tmp[:, iy:yc, ix:xc]
                    Rrs_tmp_ = ((Rcorr - Rf) / np.pi)
                    p_slope_tmp[:, iy:yc, ix:xc] = _p_slope_
                    Rrs_tmp[:, iy:yc, ix:xc] = Rrs_tmp_
                    # no sunglint removal
                    # Rrs_tmp[:, iy:yc, ix:xc] = (Rcorr  / np.pi)/ Ttot_du

            print('success')
            xres_ = xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rrs_tmp),
                                    p_slope=(['wl', "y", "x"], p_slope_tmp),
                                    BRDFg=(["y", "x"], Rf_tmp),
                                    ),
                               coords=dict(wl=raster.wl,
                                           x=raster.x,
                                           y=raster.y,
                                           aot_ref=_aot_ref),
                               ).__deepcopy__()
            xres.append(xres_)
        xres = xr.concat(xres, dim='aot_ref')
        self.xres = xres

        metric = (xres.Rrs.isel(wl=slice(0, 7)).min('wl'))
        # remove pixel with negative Rrs in Rayleigh correction (e.g., shadows)
        metric = metric.where(metric.isel(aot_ref=0) > 0)

        metric_pixnum = metric.where(metric < 0).count(['x', 'y'])
        neg_pix_max = np.max([self.neg_pix_max, metric_pixnum.max()])
        metric_pixnum_max = neg_pix_max * neg_pix_rel_thresh

        aot_ref_maxs = metric_pixnum.where(metric_pixnum < metric_pixnum_max, drop=True).aot_ref
        if len(aot_ref_maxs) > 0:
            aot_ref_max = aot_ref_maxs.values[-1]
        else:
            aot_ref_max = aot_ref_maxs
        ind = np.nanargmin((aot_refs - aot_ref_max) ** 2)

        map_pixels = metric.sel(aot_ref=aot_refs[ind])
        map_pixels = self.cams.cams_aod.interp(wl=550, method='quadratic').interp(x=metric.x, y=metric.y).where(
            map_pixels < 0)

        cams_aot_ref_mean_ = map_pixels.mean().values
        if np.isnan(cams_aot_ref_mean_):
            cams_aot_ref_mean_ = self.cams.cams_aod.interp(wl=550, method='quadratic').mean()

        self.recalib_aot_ref = float(aot_ref_max / cams_aot_ref_mean_)
        self.retrieved_aot_ref = aot_ref_max

    def final_process(self,
                      raster,
                      aot_ref_raster,
                      iwl_swir=[-2, -1],
                      chunk=512
                      ):
        ##############################################
        # Estimate AOT from visible-NIR
        ##############################################
        monoview = self.monoview
        _R_ = self._R_

        _Nwl, _height, _width = raster.bands.shape
        sunglint_eps = self.sunglint_eps.values

        sigma2 = (self.wind_img.wind + 0.586) / 195.3
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)

        pressure = self.cams.raster.sp
        Tg_diff_raster = self.Tg_diff_raster
        Tg_raster = self.Tg_raster

        Rrs_tmp = np.full((_Nwl, _height, _width), np.nan, dtype=self.prod._type)
        # Ratm = np.full((_Nwl, _height, _width), np.nan, dtype=self.prod._type)

        # p_slope_tmp = np.full((_Nwl, _height, _width), np.nan, dtype=self.prod._type)
        Rf_tmp = np.full((_height, _width), np.nan, dtype=self.prod._type)
        aot_tmp = np.full((_height, _width), np.nan, dtype=self.prod._type)
        aot_ref_bound = np.min([self.aot_ref_max * 1.2, self.aot_ref_lim])
        #aot_refs = np.linspace(0, aot_ref_bound, 21)
        Ttot_Ed_lut = self.trans_aero_lut.Ttot_Ed #.interp(aot_ref=aot_refs, method='quadratic')

        for iy in range(0, _height, chunk):
            yc = min(_height, iy + chunk)

            for ix in range(0, _width, chunk):
                xc = min(_width, ix + chunk)

                _band_rad = raster.bands[:, iy:yc, ix:xc]
                x_subset, y_subset = _band_rad.x, _band_rad.y

                Nwl, Ny, Nx = _band_rad.shape
                if Ny == 0 or Nx == 0:
                    continue
                arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

                # subsetting
                _sigma2_ = _sigma2[iy:yc, ix:xc]
                _sza = raster.sza[iy:yc, ix:xc]  # .values
                _pressure = pressure.interp(x=x_subset, y=y_subset).values
                _Tg_diff_raster = Tg_diff_raster.interp(x=x_subset, y=y_subset).values
                _Tg_raster = Tg_raster.interp(x=x_subset, y=y_subset).values
                _aot_ref_raster = aot_ref_raster.interp(x=x_subset, y=y_subset).values

                if monoview:
                    _raa = raster.raa[iy:yc, ix:xc]
                    _vza = raster.vza[iy:yc, ix:xc]
                    _vza_mean = _vza.values
                else:
                    _raa = raster.raa[:, iy:yc, ix:xc]
                    _vza = raster.vza[:, iy:yc, ix:xc]
                    _vza_mean = np.mean(_vza, axis=0).values

                _azi = (180. - _raa) % 360
                _air_mass_ = Misc.air_mass(_sza, _vza)
                _p_slope_ = self.p_slope(_sza, _vza, _raa, sigma2=_sigma2_,
                                         monoview=monoview).values

                _pressure_ = _pressure / self.pressure_ref

                # construct wl,y,x raster for Rayleigh optical thickness
                _rot_raster = _R_._multiplicate(self.rot.values, _pressure_, arr_tmp)

                # get LUT values
                _Rdiff = _R_.interp_Rlut(self.szas, _sza.values,
                                         self.vzas, _vza.values,
                                         self.azis, _azi.values,
                                         self.aot_refs, _aot_ref_raster,
                                         Nwl, Ny, Nx, self.Rdiff_lut.values)

                _Rdiff = _Rdiff * _Tg_diff_raster * _pressure_

                _aot = _R_._interp_aotlut(self.aot_lut.aot_ref.values,
                                          _aot_ref_raster, Nwl, Ny, Nx, self.aot_lut.values)

                #  correction for diffuse light
                Rcorr = _band_rad.values - _Rdiff

                # direct transmittance up/down
                Tdir = Misc.transmittance_dir(_aot, _air_mass_, _rot_raster)

                # vTotal transmittance (for Ed and Lu)
                Tdown = _R_._interp_Tlut(self.szas, _sza.values, Ttot_Ed_lut.aot_ref.values,
                                         _aot_ref_raster, Nwl, Ny, Nx,
                                         Ttot_Ed_lut.values)
                Tup = _R_._interp_Tlut(self.vzas, _vza_mean, Ttot_Ed_lut.aot_ref.values,
                                       _aot_ref_raster, Nwl, Ny, Nx, Ttot_Ed_lut.values ** 1.05)
                Ttot_du = Tdown * Tup * _Tg_raster

                self._p_slope_ = _p_slope_
                self._Tg_raster = _Tg_raster
                self._p_slope_ = _p_slope_

                Rf = np.full((len(iwl_swir), Ny, Nx), np.nan, dtype=np.float32)
                for iwl in iwl_swir:
                    if monoview:
                        Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_raster[iwl] * sunglint_eps[iwl] * _p_slope_)
                    else:
                        Rf[iwl] = (sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl])
                        Rf[iwl] = Rf[iwl] / (
                                _Tg_raster[iwl] * Tdir[iwl] * sunglint_eps[iwl] * _p_slope_[iwl])

                Rf[Rf < -0.001] = -0.00001  # np.nan
                Rf[Rf < 0.] = 0.0
                Rf = np.min(Rf, axis=0)
                Rf_tmp[iy:yc, ix:xc] = Rf

                Rf = _R_._multiplicate(sunglint_eps, Rf, arr_tmp)
                Rf = _Tg_raster * Tdir * Rf * _p_slope_ / (sunglint_eps[-1] * _p_slope_[-1])

                # sunglint removal
                # Rrs_tmp_ = Rrs_tmp[:, iy:yc, ix:xc]
                Rrs_tmp_ = ((Rcorr - Rf) / np.pi)
                # p_slope_tmp[:, iy:yc, ix:xc] = _p_slope_

                # Convert from TOA to BOA for positive values
                Ttot_du[Rrs_tmp_ < 0] = 1.
                Rrs_tmp_ = Rrs_tmp_ / Ttot_du
                Rrs_tmp[:, iy:yc, ix:xc] = Rrs_tmp_
                aot_tmp[iy:yc, ix:xc] = _aot_ref_raster

                # no sunglint removal
                # Rrs_tmp[:, iy:yc, ix:xc] = (Rcorr  / np.pi)/ Ttot_du

        xres = xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rrs_tmp),
                               # p_slope=(['wl', "y", "x"], p_slope_tmp),
                               BRDFg=(["y", "x"], Rf_tmp),
                               aot550=(["y", "x"], aot_tmp)
                               ),
                          coords=dict(wl=raster.wl,
                                      x=raster.x,
                                      y=raster.y,
                                      ),
                          )

        self.xres = xres

    @staticmethod
    @njit
    def filter2d(image, weight, windows):
        '''
         Function to convolve parameter image with uncertainty image
        :param image: parameter image
        :param weight: uncertainty image
        :param windows: size of the window for convolution
        :return: convolved result with same shape as image

        '''
        M, N = np.shape(image)
        Mf, Nf = windows
        Mf2 = Mf // 2
        Nf2 = Nf // 2
        threshold = 0
        result = image
        for i in range(M):
            for j in range(N):
                num = 0.0
                norm = 0.0
                if weight[i, j] > threshold:
                    for ii in range(Mf):
                        ix = i - Mf2 + ii
                        if ix < M:
                            for jj in range(Nf):

                                iy = j - Nf2 + jj
                                if iy < N:
                                    wgt = weight[ix, iy]
                                    if wgt > 0.:
                                        num += (wgt * image[ix, iy])
                                        norm += wgt
                    result[i, j] = num / norm
        return result

    @staticmethod
    def conv_mapping(x):
        """
        Nan-mean convolution
        """
        # get index of central pixel
        idx = len(x) // 2
        if np.isnan(x[idx]) and not np.isnan(np.delete(x, idx)).all():
            return np.nanmean(np.delete(x, idx))
        elif np.isnan(np.delete(x, idx)).all():
            return x[idx]
        else:
            return np.nanmean(x)

    def smoothing(self,
                  raster,
                  varname='wind',
                  windows=np.array([15, 15]),
                  mask=np.ones((5, 5))):
        weights = 1 + 0 * (1 / raster[varname] ** 2).values
        param = raster[varname].values
        raster_smoothed = self.filter2d(param, weights, windows)
        res = ndimage.generic_filter(raster_smoothed, function=self.conv_mapping, footprint=mask, mode='nearest')

        return xr.Dataset({varname: (["y", "x"], res)}, coords=dict(y=raster.y, x=raster.x))
