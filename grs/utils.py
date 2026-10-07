'''
Atmospheric correction utilities to manage LUTs and atmosphere parameters (aerosols, gases).

Contents
========

``Rasterization``
    Numba kernels that fill a (wl, y, x) raster from a look-up table (LUT) by nearest-neighbour
    lookup on per-pixel geometry (solar/view zenith, relative azimuth) and aerosol optical thickness.
``Aerosol``
    Aerosol parameters, spectral AOT fit and fine/coarse mode mixing fit.
``Gases`` / ``GaseousTransmittance``
    Absorbing gas parameters, and band-averaged direct gaseous transmittance rasters computed
    from hyperspectral absorption LUTs and CAMS total columns.
``Misc``
    Small helpers (pressure vs altitude, direct transmittance, air mass, downsampling).

Performance notes
=================

* The nearest-index search no longer allocates an array per call, the kernels run in parallel
  (``prange``) and are cached on disk (``cache=True``): the first call in a new environment compiles,
  later runs reuse the cache. Index maps are computed once per pixel, not once per band.
* The band-averaged transmittance is computed as a weighted sum with pre-computed spectral-response
  weights (one ``exp`` table and one dot product per band) instead of a Python loop over
  concentrations x bands with xarray selections and trapezoid integrals.
* The spectral AOT fit and the fine-mode fit are linear least-squares problems and are solved in
  closed form (``numpy``), so ``scipy`` is no longer needed.
* The parallel kernels must be called from a single thread at a time (numba's default threading layer
  is not safe for concurrent parallel calls); do not call them from several dask threads at once.
'''
import logging

import numpy as np
import xarray as xr
import rioxarray

from numba import njit, prange


@njit(cache=True)
def _getnearpos(array, value):
    '''
    Get the index of the nearest neighbour of `value` in a 1-D array.

    Allocation-free linear scan (the previous ``np.abs(array - value).argmin()`` created a temporary
    array at every call). Same result, including ties (first index) and NaN `value` (index 0).

    :param array: 1-D array to search in (need not be sorted)
    :param value: value whose nearest index is wanted
    :return: index of the element of `array` closest to `value`
    '''
    best = 0
    bestd = abs(array[0] - value)
    for i in range(1, array.shape[0]):
        d = abs(array[i] - value)
        if d < bestd:
            bestd = d
            best = i
    return best


@njit(parallel=True, cache=True)
def _nearest_map(axis, values):
    '''
    Nearest-neighbour index of every pixel of a 2-D array on a 1-D axis.

    :param axis: 1-D LUT axis
    :param values: 2-D array (Ny, Nx) of values
    :return: int32 array (Ny, Nx) of indices into `axis`, -1 where `values` is NaN
    '''
    Ny, Nx = values.shape
    out = np.full((Ny, Nx), -1, dtype=np.int32)
    for iy in prange(Ny):
        for ix in range(Nx):
            v = values[iy, ix]
            if not np.isnan(v):
                out[iy, ix] = _getnearpos(axis, v)
    return out



class Rasterization:
    '''
    Class holding numba-compiled kernels that sample a LUT onto an image grid.

    All kernels take the LUT axes (``szas``, ``vzas``, ``azis``, ``aot_refs``) and the per-pixel values
    of the corresponding quantity (``_sza``, ``_vza``, ``_azi``, ``_aot_ref``), pick the nearest LUT node
    and return a float32 array ``(Nwl, Ny, Nx)``. A pixel is NaN in the output as soon as one of its
    inputs is NaN (before, only a NaN solar zenith angle was tested and other NaN inputs silently
    selected node 0).

    LUT shapes: Rayleigh LUT ``(Nwl, Nsza, Nvza, Nazi)``, aerosol LUT ``(Naot, Nwl, Nsza, Nvza, Nazi)``,
    transmittance LUT ``(Naot, Nwl, Nsza)``, AOT LUT ``(Naot, Nwl)``.
    '''

    def __init__(self,
                 monoview=False):
        '''
        Set compiled function for:
            - monoview = True : same viewing angles for all the spectral bands
              (``_vza`` and ``_azi`` are 2-D arrays (Ny, Nx))
            - monoview = False : viewing angles depend on spectral band, e.g. Sentinel-2 images
              (``_vza`` and ``_azi`` are 3-D arrays (Nwl, Ny, Nx))

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
    @njit(parallel=True, cache=True)
    def _interp_Rlut_rayleigh_mono(szas, _sza,
                                   vzas, _vza,
                                   azis, _azi,
                                   Nwl, Ny, Nx, lut):
        '''Rayleigh reflectance LUT ``lut[Nwl, Nsza, Nvza, Nazi]``, same geometry for all bands.'''
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        isza = _nearest_map(szas, _sza)
        ivza = _nearest_map(vzas, _vza)
        iazi = _nearest_map(azis, _azi)
        for _iwl in prange(Nwl):
            for _iy in range(Ny):
                for _ix in range(Nx):
                    a = isza[_iy, _ix]
                    b = ivza[_iy, _ix]
                    c = iazi[_iy, _ix]
                    if a < 0 or b < 0 or c < 0:
                        continue
                    arr_lut[_iwl, _iy, _ix] = lut[_iwl, a, b, c]
        return arr_lut

    @staticmethod
    @njit(parallel=True, cache=True)
    def _interp_Rlut_mono(szas, _sza,
                          vzas, _vza,
                          azis, _azi,
                          aot_refs, _aot_ref,
                          Nwl, Ny, Nx, lut):
        '''Aerosol reflectance LUT ``lut[Naot, Nwl, Nsza, Nvza, Nazi]``, same geometry for all bands.'''
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        isza = _nearest_map(szas, _sza)
        ivza = _nearest_map(vzas, _vza)
        iazi = _nearest_map(azis, _azi)
        iaot = _nearest_map(aot_refs, _aot_ref)
        for _iwl in prange(Nwl):
            for _iy in range(Ny):
                for _ix in range(Nx):
                    a = isza[_iy, _ix]
                    b = ivza[_iy, _ix]
                    c = iazi[_iy, _ix]
                    d = iaot[_iy, _ix]
                    if a < 0 or b < 0 or c < 0 or d < 0:
                        continue
                    arr_lut[_iwl, _iy, _ix] = lut[d, _iwl, a, b, c]
        return arr_lut

    @staticmethod
    @njit(parallel=True, cache=True)
    def _interp_Rlut_rayleigh(szas, _sza,
                              vzas, _vza,
                              azis, _azi,
                              Nwl, Ny, Nx, lut):
        '''Rayleigh reflectance LUT ``lut[Nwl, Nsza, Nvza, Nazi]``, band-dependent view angles.'''
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        isza = _nearest_map(szas, _sza)
        for _iwl in prange(Nwl):
            for _iy in range(Ny):
                for _ix in range(Nx):
                    a = isza[_iy, _ix]
                    vz = _vza[_iwl, _iy, _ix]
                    az = _azi[_iwl, _iy, _ix]
                    if a < 0 or np.isnan(vz) or np.isnan(az):
                        continue
                    ivza = _getnearpos(vzas, vz)
                    iazi = _getnearpos(azis, az)
                    arr_lut[_iwl, _iy, _ix] = lut[_iwl, a, ivza, iazi]
        return arr_lut

    @staticmethod
    @njit(parallel=True, cache=True)
    def _interp_Rlut(szas, _sza,
                     vzas, _vza,
                     azis, _azi,
                     aot_refs, _aot_ref,
                     Nwl, Ny, Nx, lut):
        '''Aerosol reflectance LUT ``lut[Naot, Nwl, Nsza, Nvza, Nazi]``, band-dependent view angles.'''
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        isza = _nearest_map(szas, _sza)
        iaot = _nearest_map(aot_refs, _aot_ref)
        for _iwl in prange(Nwl):
            for _iy in range(Ny):
                for _ix in range(Nx):
                    a = isza[_iy, _ix]
                    d = iaot[_iy, _ix]
                    vz = _vza[_iwl, _iy, _ix]
                    az = _azi[_iwl, _iy, _ix]
                    if a < 0 or d < 0 or np.isnan(vz) or np.isnan(az):
                        continue
                    ivza = _getnearpos(vzas, vz)
                    iazi = _getnearpos(azis, az)
                    arr_lut[_iwl, _iy, _ix] = lut[d, _iwl, a, ivza, iazi]
        return arr_lut

    @staticmethod
    @njit(parallel=True, cache=True)
    def _interp_Tlut(szas, _sza,
                     aot_refs, _aot_ref,
                     Nwl, Ny, Nx, lut):
        '''Transmittance LUT ``lut[Naot, Nwl, Nsza]``.'''
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        isza = _nearest_map(szas, _sza)
        iaot = _nearest_map(aot_refs, _aot_ref)
        for _iwl in prange(Nwl):
            for _iy in range(Ny):
                for _ix in range(Nx):
                    a = isza[_iy, _ix]
                    d = iaot[_iy, _ix]
                    if a < 0 or d < 0:
                        continue
                    arr_lut[_iwl, _iy, _ix] = lut[d, _iwl, a]
        return arr_lut

    @staticmethod
    @njit(parallel=True, cache=True)
    def _interp_aotlut(aot_refs, _aot_ref,
                       Nwl, Ny, Nx, lut):
        '''AOT LUT ``lut[Naot, Nwl]`` (spectral scaling of the reference AOT).'''
        arr_lut = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)
        iaot = _nearest_map(aot_refs, _aot_ref)
        for _iwl in prange(Nwl):
            for _iy in range(Ny):
                for _ix in range(Nx):
                    d = iaot[_iy, _ix]
                    if d < 0:
                        continue
                    arr_lut[_iwl, _iy, _ix] = lut[d, _iwl]
        return arr_lut

    @staticmethod
    @njit(parallel=True, cache=True)
    def _multiplicate(arrwl, raster, arresult):
        '''
                Fill ``arresult[iwl] = arrwl[iwl] * raster`` for every band, in place and without temporaries.

                :param arrwl: per-band factor, shape (Nwl,) (or (Nwl, Ny, Nx))
                :param raster: 2-D array (Ny, Nx)
                :param arresult: output array (Nwl, Ny, Nx), overwritten and returned
                '''
        Nwl, Ny, Nx = arresult.shape
        for iwl in prange(Nwl):
            np.multiply(arrwl[iwl], raster, arresult[iwl])
        return arresult


class Aerosol:
    '''
    Aerosol parameters and parameterizations.

    Attributes
    ----------
    aot550 : float
        Aerosol optical thickness at 550 nm (default 0.1).
    wavelengths, aot, wl, ang : list
        Containers for measured/derived spectral AOT data and Angstrom exponents (set by the user).
    o3du, no2du : float
        Ozone and NO2 total columns in Dobson units (defaults 300).
    fcoef : float or numpy.ndarray
        Fraction of the fine mode in the bimodal aerosol model (0.5 by default;
        a 1-element array after `fit_aero`).
    popt : numpy.ndarray
        Coefficients ``[a, b, c]`` of the spectral AOT fit (empty until `fit_spectral_aot` is called).
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
        '''Second-order polynomial in ln(wavelength): spectral variation of ln(AOT).'''

        return (a + b * lnwl + c * lnwl ** 2)

    def fit_spectral_aot(self, wl, aot):
        '''
        Fit ``ln(aot) = a + b ln(wl) + c ln(wl)^2`` and store ``[a, b, c]`` in `popt`.

        The model is linear in its coefficients, so the least-squares solution is obtained directly
        with ``numpy.polyfit`` (same result as an iterative ``curve_fit``).

        :param wl: wavelengths (at least 3)
        :param aot: aerosol optical thickness at `wl` (all > 0)
        '''
        wl = np.asarray(wl, dtype=float)
        aot = np.asarray(aot, dtype=float)
        if wl.size < 3:
            raise ValueError('at least 3 wavelengths are needed to fit a second-order spectral AOT model')
        if np.any(aot <= 0) or np.any(wl <= 0):
            raise ValueError('wavelengths and AOT must be strictly positive')
        c, b, a = np.polyfit(np.log(wl), np.log(aot), 2)
        self.popt = np.array([a, b, c])

    def get_spectral_aot(self, wl):
        '''
        AOT at the wavelengths `wl` from the fitted spectral model.

        :param wl: wavelengths
        :return: AOT array (requires a previous call to `fit_spectral_aot`)
        '''
        if len(self.popt) == 0:
            raise RuntimeError('call fit_spectral_aot first')
        lnwl = np.log(wl)
        return np.exp(self.func(lnwl, *self.popt))

    def func_aero(self, Cext, fcoef):
        '''
        Bimodal aerosol model: normalized spectral extinction as a mix of fine and coarse modes,
        ``fcoef * Cext[0] + (1 - fcoef) * Cext[1]``.
        '''
        return fcoef * Cext[0] + (1 - fcoef) * Cext[1]

    def fit_aero(self, nCext_f, nCext_c, naot):
        '''
        Fine-mode fraction from a fit onto the normalized spectral AOT (e.g. from AERONET).

        The model is linear in the fraction, so the least-squares solution is exact and computed
        directly: ``f = sum(d * (naot - Cc)) / sum(d * d)`` with ``d = Cf - Cc``.

        Parameters
        ----------
        nCext_f : array-like
            Normalized extinction coefficient of the fine mode aerosols.
        nCext_c : array-like
            Normalized extinction coefficient of the coarse mode aerosols.
        naot : array-like
            Normalized spectral aerosol optical thickness.

        Returns
        -------
        numpy.ndarray
            One-element array with the fine mode fraction (also stored in `fcoef`). The value is not
            constrained to [0, 1].
        '''
        cf = np.asarray(nCext_f, dtype=float)
        cc = np.asarray(nCext_c, dtype=float)
        d = cf - cc
        self.fcoef = np.array([np.sum(d * (np.asarray(naot, dtype=float) - cc)) / np.sum(d * d)])
        return self.fcoef


class CamsParams:
    '''
        Description of a CAMS variable.

        :param name: variable name in the CAMS dataset (e.g. ``'tcwv'``)
        :param resol: number of decimals to which the values are rounded (limits the number of distinct
            concentrations for which a transmittance is computed)
    '''

    def __init__(self,
                 name,
                 resol):
        self.name = name
        self.resol = resol


class Gases():
    '''
    Intermediate class to set parameters for absorbing gases.

    Attributes
    ----------
    pressure : float
        Surface pressure (hPa).
    psl : float
        Sea level pressure (hPa).
    to3c, tno2c, tch4c : float
        Default total columns of O3, NO2 and CH4 used by the obsolete `get_gaseous_optical_thickness`
        (TODO in the code: get them from CAMS).
    coef_abs_scat : dict
        Per-gas factor (6S-like) accounting for the vertical distribution of the absorber relative to a
        scattering layer above it. It scales the CAMS column of ch4, no2, o3, h2o in `Tgas`, and the
        co2, o2, o4 absorption of the background mixture in `Tgas_background`.
    '''

    def __init__(self):
        # atmosphere auxiliary data
        # TODO get them from CAMS
        self.pressure = 1010
        self.to3c = 6.5e-3
        self.tno2c = 3e-6
        self.tch4c = 1e-2
        self.psl = 1013
        self.coef_abs_scat = {'co2':0.4,
                              'o2':0.3,
                              'o4':0.3,
                              'ch4': 0.5,
                              'no2': 1,
                              'o3': 1,
                              'h2o':0.3}

def _srf_supports(srf, wl):
    '''
    Integration weights of the spectral response functions (SRF).

    For each band the weights ``w`` satisfy ``sum(w * T) == trapz(T * srf, wl) / trapz(srf, wl)``
    where the integrals run over the wavelengths at which the SRF is defined (not NaN), exactly as
    a trapezoid integration of the selected samples.

    :param srf: array (Nbands, Nwl), NaN outside the support of each band
    :param wl: 1-D array (Nwl) of wavelengths (monotonic)
    :return: list of (indices, weights) per band; weights are NaN if the band has less than 2 samples
    '''
    supports = []
    for s in np.asarray(srf, dtype=float):
        idx = np.flatnonzero(~np.isnan(s))
        if idx.size < 2:
            supports.append((idx, np.full(idx.size, np.nan)))
            continue
        x = wl[idx]
        dx = np.diff(x)
        c = np.zeros(idx.size)
        c[:-1] += dx / 2
        c[1:] += dx / 2
        w = s[idx] * c
        supports.append((idx, w / w.sum()))
    return supports


def _band_average(k, air_mass, vals, supports, max_elements=20_000_000):
    '''
    Band-averaged direct transmittance ``exp(-air_mass * k * val)`` for several amounts `val`.

    :param k: 1-D hyperspectral absorption coefficient (Nwl)
    :param air_mass: scalar air mass
    :param vals: 1-D array (Nvals) of gas amounts (or pressures)
    :param supports: output of `_srf_supports`
    :param max_elements: size limit of the temporary exp table (the values are processed in blocks)
    :return: array (Nbands, Nvals)
    '''
    vals = np.asarray(vals, dtype=float)
    out = np.empty((len(supports), vals.size))
    for b, (idx, w) in enumerate(supports):
        kb = -air_mass * k[idx]
        step = max(1, int(max_elements // max(idx.size, 1)))
        for v0 in range(0, vals.size, step):
            out[b, v0:v0 + step] = w @ np.exp(np.outer(kb, vals[v0:v0 + step]))
    return out


class GaseousTransmittance(Gases):
    '''
    Class containing functions to compute rasters of the direct transmittance of the absorbing gases.

    Parameters
    ----------
    prod : object
        Product with ``raster`` (dataset containing the spectral response functions ``SRF`` with
        dimensions (wl, wl_hr)), ``gas_lut`` (hyperspectral absorption coefficients of co, co2, o2, o4,
        ch4, no2, o3, h2o on ``wl``), ``Twv_lut`` and the scalar ``air_mass_mean``.
    cams : object
        CAMS product with ``raster`` containing ``sp`` (surface pressure, Pa), ``tc_ch4``, ``tcno2``,
        ``gtco3`` and ``tcwv``.

    Notes
    -----
    The transmittance of a band is the SRF-weighted mean of the hyperspectral transmittance. Because
    CAMS values are rounded, it is computed once for each distinct value (a small table) and then
    interpolated (exactly, at table nodes) onto the CAMS grid.
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
        self.pressure = cams.raster.sp * 1e-2 # Pa -> hPa
        #self.coef_abs_scat = 0.3
        self.Tg_tot_coarse = None
        self.Twv_raster = None
        self._srf_cache = None
        self.cams_gases = {'ch4': CamsParams('tc_ch4', 4),
                           'no2': CamsParams('tcno2', 7),
                           'o3': CamsParams('gtco3', 4),
                           'h2o': CamsParams('tcwv', 1), }

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _srf_data(self):
        '''
        SRFs interpolated on the LUT wavelengths, as integration weights (computed once and cached).

        :return: (supports, wl_bands): per-band (indices, weights) and the band wavelengths
        '''
        if self._srf_cache is None:
            wl_ref = self.gas_lut.wl
            SRF_hr = self.prod.raster.SRF.interp(wl_hr=wl_ref.values)
            srf = SRF_hr.transpose('wl', 'wl_hr').values
            supports = _srf_supports(srf, SRF_hr.wl_hr.values)
            self._srf_cache = (supports, SRF_hr.wl.values)
        return self._srf_cache

    @staticmethod
    def _unique_values(array):
        '''Sorted distinct non-NaN values of `array`; a single value is paired with 1.2 x itself
        (the later interpolation needs at least two nodes).'''
        vals = np.unique(np.asarray(array))
        vals = vals[~np.isnan(vals)]
        if vals.size == 0:
            raise ValueError('no valid value to compute the transmittance for')
        if vals.size == 1:
            vals = np.unique(np.concatenate([vals, 1.2 * vals]))
        return vals

    def _band_table(self, k, vals, dim):
        '''
        Table of band-averaged transmittance for the amounts `vals`.

        :param k: hyperspectral absorption coefficient (on the LUT wavelengths)
        :param vals: 1-D array of amounts
        :param dim: name of the amount dimension of the result (e.g. ``'tc'``)
        :return: DataArray ``Ttot`` with dimensions (dim, wl)
        '''
        supports, wl_bands = self._srf_data()
        table = _band_average(np.asarray(k, dtype=float), float(self.air_mass_mean), vals, supports)
        return xr.DataArray(table.T, dims=(dim, 'wl'),
                            coords={dim: vals, 'wl': wl_bands}, name='Ttot')

    # ------------------------------------------------------------------
    # transmittances
    # ------------------------------------------------------------------
    def Tgas_background(self):
        '''
        Compute direct transmittance for the background absorbing gases: :math:`CO,\\ CO_2,\\ O_2,\\ O_4`,
        for the surface pressure of the CAMS grid (rounded to 0.1 hPa).

        :return: DataArray of the band transmittance with the dimensions of the pressure raster plus ``wl``
        '''
        gl = self.gas_lut
        pressure = self.pressure.round(1)
        self.ot_air = (gl.co + self.coef_abs_scat['co2'] * gl.co2 +
                       self.coef_abs_scat['o2'] * gl.o2 +
                       self.coef_abs_scat['o4'] * gl.o4) / 1000
        vals = self._unique_values(pressure)
        Tg_raster = self._band_table(self.ot_air.values, vals, 'pressure')
        return Tg_raster.interp(pressure=pressure).drop_vars(['pressure'])

        # wl_ref = gl.wl
        # SRF_hr = self.prod.raster.SRF.interp(wl_hr=wl_ref.values)
        # vals = np.unique(pressure)
        # vals = vals[~np.isnan(vals)]
        # if len(vals) == 1:
        #     vals = np.concatenate([vals, 1.2 * vals])
        # Tg_raster = []
        # for val in vals:
        #     Tg = np.exp(- self.air_mass_mean * self.ot_air * val)
        #     Tg = Tg.rename({'wl': 'wl_hr'})
        #
        #     Tg_int = []
        #     for label, srf in SRF_hr.groupby('wl', squeeze=False):
        #         srf = srf.dropna('wl_hr').squeeze()
        #         Tg_ = Tg.sel(wl_hr=srf.wl_hr)
        #         wl_integr = Tg_.wl_hr.values
        #
        #         Tg_ = np.trapezoid(Tg_ * srf, wl_integr) / np.trapezoid(srf, wl_integr)
        #         Tg_int.append(Tg_)
        #     Tg_raster.append(xr.DataArray(Tg_int, name='Ttot', coords={'wl': SRF_hr.wl.values}
        #                                   ).assign_coords({'pressure': val}))
        # Tg_raster = xr.concat(Tg_raster, dim='pressure')
        # return Tg_raster.interp(pressure=pressure).drop_vars(['pressure'])

    def Tgas(self,
             gas_name,
             coef_abs_scat=1):
        '''
        Compute hyperspectral transmittance for a given absorbing gas and
        convolve it with the spectral response functions of the satellite sensor.

        :param gas_name: name of the absorbing gas, one of ``'ch4'``, ``'no2'``, ``'o3'``, ``'h2o'``
        :param coef_abs_scat: factor applied to the CAMS column to account for the vertical distribution
            of the gas relative to the scattering layer (6S approach)
        :return: gaseous transmittance for the satellite bands, on the CAMS grid
        '''


        cams_gas = self.cams_gases[gas_name].name
        resol = self.cams_gases[gas_name].resol
        lut_abs = self.gas_lut[gas_name]

        # round number to speed up computation
        # and scale the concentration following 6S approach due to vertical distribution
        # with scattering layer above absorbing gases
        rounded = coef_abs_scat * self.cams.raster[cams_gas].round(resol)

        vals = self._unique_values(rounded)
        Tg_raster = self._band_table(lut_abs.values, vals, 'tc')
        return Tg_raster.interp(tc=rounded)

    def get_gaseous_optical_thickness(self):
        '''
        Get gaseous optical thickness from total column integrated concentration.
        Obsolete: only used by `get_gaseous_transmittance_old`. Sets ``abs_gas_opt_thick``.
        '''

        gas_lut = self.gas_lut

        ot_o3 = gas_lut.o3 * self.to3c
        ot_ch4 = gas_lut.ch4 * self.tch4c
        ot_no2 = gas_lut.no2 * self.tno2c
        ot_air = (gas_lut.co + self.coef_abs_scat['co2'] * gas_lut.co2 +
                  self.coef_abs_scat['o2'] * gas_lut.o2 +
                  self.coef_abs_scat['o4'] * gas_lut.o4) * self.pressure / 1000
        self.abs_gas_opt_thick = ot_ch4 + ot_no2 + ot_o3 + ot_air

    def get_gaseous_transmittance(self,
                                  gases=['ch4','no2','o3','h2o'],
                                  background=True):
        '''
        Get the final total gaseous transmittance: product of the transmittances of `gases`
        (CAMS columns) and, optionally, of the background gases.

        :param gases: gases to include, among ``'ch4'``, ``'no2'``, ``'o3'``, ``'h2o'``
        :param background: also include `Tgas_background`
        :return: DataArray of the total transmittance on the (coarse) CAMS grid; also stored in
            ``Tg_tot_coarse``
        '''


        first = True
        for gas in gases:
            if first:
                Tg_tot = self.Tgas(gas, self.coef_abs_scat[gas])
                first=False
            else:
                Tg_tot =Tg_tot *self.Tgas(gas,
                                      self.coef_abs_scat[gas])

        if background:
            Tg_tot = Tg_tot * self.Tgas_background()

        self.Tg_tot_coarse = Tg_tot
        # TODO remove interp for the whole object and proceed with loop on spectral bands to save memory
        return Tg_tot  # .interp(x=self.prod.raster.x, y=self.prod.raster.y)

    def correct_gaseous_transmittance(self):

        return

    def get_gaseous_transmittance_old(self):
        '''
        Obsolete function: band transmittance of all gases together from the default total columns
        (no CAMS data). Sets ``Tg_other``.
        '''
        self.get_gaseous_optical_thickness()
        supports, wl_bands = self._srf_data()
        # exp(-air_mass * sum(k_i * c_i)) = _band_average with val = 1 and k = total optical thickness
        k = np.asarray(self.abs_gas_opt_thick.values, dtype=float)
        table = _band_average(k, float(self.air_mass_mean), np.array([1.0]), supports)
        self.Tg_other = xr.DataArray(table[:, 0], name='Ttot', coords={'wl': wl_bands})

    def other_gas_correction(self, raster_name='masked_raster', variable='Rtoa_masked'):
        '''
        Correct for transmittance of high altitude gases
        (i.e., no coupling with low atmosphere scattering).

        Divides ``self.<raster_name>[variable]`` by ``Tg_other`` (computed with
        `get_gaseous_transmittance_old` if needed) and flags the variable with the attribute
        ``other_gas_correction``.

        :param raster_name: name of the attribute holding the raster dataset
        :param variable: variable to correct
        '''
        raster = self.__dict__[raster_name]
        attrs = raster[variable].attrs
        if attrs.__contains__('other_gas_correction'):
            if attrs['other_gas_correction']:
                print('raster ' + raster_name + '.' + variable + ' is already corrected for other gases transmittance')
                print('set attribute other_gas_correction to False to proceed anyway')
                return
        if self.Tg_other is None:
            # TODO check and replace obsolete function
            self.get_gaseous_transmittance_old(self.air_mass_mean)
        raster[variable] = raster[variable] / self.Tg_other
        raster[variable].attrs['other_gas_correction'] = True

    def water_vapor_correction(self, raster_name='coarse_masked_raster', variable='Rtoa_masked'):
        '''
        Correct for water vapor transmittance: divides ``self.<raster_name>[variable]`` by the
        water vapor transmittance raster (see `get_wv_transmittance_raster`) and flags the variable with
        the attribute ``water_vapor_correction``.

        :param raster_name: name of the attribute holding the raster dataset
        :param variable: variable to correct
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
        Get the water vapor transmittance raster for correction of the image (stored in ``Twv_raster``).

        :param tcwv_raster: dataset with the total column water vapor ``tcwv``
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
        '''Compute the pressure for a given altitude (barometric formula; NaN altitudes are
        treated as sea level).

        :param alt: altitude in meters (float or np.array)
        :param psl: pressure at sea level in hPa
        :return: pressure at the given altitude in hPa'''

        palt = psl * (1. - 0.0065 * np.nan_to_num(alt) / 288.15) ** 5.255
        return palt

    @staticmethod
    def transmittance_dir(aot, air_mass, rot=0):
        r'''Direct transmittance ``exp(-(rot + aot) * air_mass)``:

        .. math::

           T_{dir} = \exp\left[-(\tau_R + \tau_a)\, M\right]

        :param aot: aerosol optical thickness
        :param air_mass: air mass
        :param rot: Rayleigh optical thickness (default 0)'''
        return np.exp(-(rot + aot) * air_mass)

    @staticmethod
    def air_mass(sza, vza):
        r'''Two-way (sun-surface-sensor) geometric air mass ``1/cos(vza) + 1/cos(sza)``; angles in degrees:

        .. math::

           M = \frac{1}{\cos\theta_v} + \frac{1}{\cos\theta_s}
        '''
        return 1 / np.cos(np.radians(vza)) + 1 / np.cos(np.radians(sza))

    @staticmethod
    def _downsample(ds, dims=('y', 'x'), factor=2, categorical=('flags', 'mask')):
        """Reduce the spatial resolution of `ds` by an integer `factor`.

        Continuous variables are block-averaged ignoring NaN. Categorical
        variables (``flags``, ``mask``) are taken from the top-left pixel of
        each block (nearest neighbour), which preserves their values and
        dtype. Edge pixels that do not fill a complete block are trimmed, so
        the origin is unchanged and the pixel size is multiplied by `factor`.

        Parameters
        ----------
        ds : xarray.Dataset
            Dataset with the spatial dimensions `dims`.
        dims : tuple of str, optional
            Names of the (y, x) dimensions.
        factor : int, optional
            Block size (default 2).
        categorical : tuple of str, optional
            Variables resampled by nearest neighbour.

        Returns
        -------
        xarray.Dataset
        """
        window = {d: factor for d in dims}
        continuous = [v for v in ds.data_vars if v not in categorical]
        coarse = ds[continuous].coarsen(window, boundary='trim').mean()
        for v in continuous:
            coarse[v].attrs = ds[v].attrs

        for v in categorical:
            if v in ds.data_vars:
                sub = ds[v].isel({d: slice(0, coarse.sizes[d] * factor, factor) for d in dims})
                sub = sub.assign_coords({d: coarse[d] for d in dims})
                coarse[v] = sub
                coarse[v].attrs = ds[v].attrs

        coarse.attrs = ds.attrs
        return coarse

    @staticmethod
    def _update_geotransform(ds, dims=('y', 'x')):
        """Rewrite the ``GeoTransform`` of the ``spatial_ref`` variable from the coordinates.

        Needed after resampling so that GIS readers see the right pixel size.
        Failures (e.g. no CRS or non-georeferenced data) only log a warning.
        """
        try:
            ds = ds.rio.set_spatial_dims(x_dim=dims[1], y_dim=dims[0])
            ds = ds.rio.write_transform(ds.rio.transform(recalc=True))
        except Exception as e:
            logging.warning(f'could not update geotransform: {e}')
        return ds

