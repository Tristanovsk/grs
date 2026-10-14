r'''
Atmospheric correction kernel: wind / aerosol retrieval and remote-sensing reflectance (dask-friendly).

Physical model
==============

Notation
--------

.. list-table::
   :widths: 25 75

   * - :math:`\lambda`
     - central wavelength of the band (``wl_true``)
   * - :math:`\theta_s,\ \theta_v,\ \phi`
     - solar zenith, viewing zenith and relative azimuth angles (``sza``, ``vza``, ``raa``); the LUTs are
       read with :math:`\phi_{LUT} = (180^\circ - \phi) \bmod 360^\circ`
   * - :math:`R_{TOA}(\lambda)`
     - top-of-atmosphere reflectance of the L1C image (``bands``)
   * - :math:`\tau_{ref}`
     - aerosol optical thickness (AOT) at the reference wavelength (550 nm, ``aot_ref``)
   * - :math:`\tau_a(\lambda),\ \tau_R(\lambda)`
     - aerosol and Rayleigh optical thicknesses (``aot_lut``, ``rot``)
   * - :math:`P,\ P_0`
     - CAMS surface pressure and reference pressure :math:`P_0 = 101320` Pa
   * - :math:`T_g,\ T_g^{diff}`
     - gaseous transmittance (H\ :sub:`2`\ O, CH\ :sub:`4`) of the direct and of the diffuse light
       (see `Kernel.set_gas_transmittance`)
   * - :math:`w_m`
     - weight of the aerosol model :math:`m` in the mixture (``aerosol_combination``)
   * - :math:`U`
     - wind speed (m s\ :sup:`-1`)
   * - :math:`\varepsilon(\lambda)`
     - mean spectral shape of the sunglint reflectance at small angles (``sunglint_eps``)
   * - :math:`P_{slope}(\lambda)`
     - probability of the wave slopes reflecting the sun towards the sensor (`Kernel.p_slope`)

Atmospheric diffuse light
-------------------------

The LUT intensities :math:`I_m` of each aerosol model are mixed and converted into reflectance
(`Kernel.lut_preparation`):

.. math::

   R_{diff}(\lambda; \tau_{ref}, \theta_s, \theta_v, \phi) =
   \frac{1}{\cos\theta_s} \sum_m w_m\, I_m(\lambda; \tau_{ref}, \theta_s, \theta_v, \phi)

and removed from the TOA signal, after scaling for the surface pressure and the gaseous absorption:

.. math::

   R_{corr}(\lambda) = R_{TOA}(\lambda) - T_g^{diff}(\lambda)\, \frac{P}{P_0}\,
   R_{diff}(\lambda; \tau_{ref}, \theta_s, \theta_v, \phi)

Direct and total transmittances
-------------------------------

With the two-way geometric air mass :math:`M = 1/\cos\theta_s + 1/\cos\theta_v`, the direct
(sun-surface-sensor) transmittance is

.. math::

   T_{dir}(\lambda) = \exp\left[-\left(\tau_a(\lambda) + \tau_R(\lambda)\,\frac{P}{P_0}\right) M\right],
   \qquad \tau_a(\lambda) = \sum_m w_m\, \tau_{a,m}(\lambda; \tau_{ref})

and the total (direct + diffuse) transmittance used to bring the water signal to the surface is

.. math::

   T_{tot}(\lambda) = T_{\downarrow}(\lambda; \theta_s, \tau_{ref})\;
   T_{\downarrow}(\lambda; \bar\theta_v, \tau_{ref})^{1.05}\; T_g(\lambda)

where :math:`T_{\downarrow}` is the LUT transmittance of the downward irradiance (``Ttot_Ed``) and
:math:`\bar\theta_v` the viewing angle averaged over the bands.

Sunglint
--------

The sea surface is described by the isotropic Cox-Munk distribution of slopes, whose mean square slope
depends on the wind speed (Cox and Munk, 1954):

.. math::

   \sigma^2 = 0.003 + 0.00512\, U = \frac{U + 0.586}{195.3}

The sunglint is estimated in a SWIR band :math:`\lambda_s`, where the water-leaving signal is null, and
extrapolated to all bands with the spectral shape :math:`\varepsilon(\lambda)` and the slope probability of
each band (the viewing angles of Sentinel-2 depend on the band):

.. math::

   R_f = \frac{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)\, R_{corr}(\lambda_s)}
              {T_g(\lambda_s)\, T_{dir}(\lambda_s)\, \varepsilon(\lambda_s)\, P_{slope}(\lambda_s)},
   \qquad
   R_g(\lambda) = T_g(\lambda)\, T_{dir}(\lambda)\, \varepsilon(\lambda)\, P_{slope}(\lambda)\,
   \frac{R_f}{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)}

with :math:`\lambda_r` a reference SWIR band (see `Kernel.aerosol_visible_chunk` and
`Kernel.final_process_chunk` for the bands used).

Remote-sensing reflectance
--------------------------

.. math::

   R_{rs}(\lambda) = \frac{R_{corr}(\lambda) - R_g(\lambda)}{\pi\, T_{tot}(\lambda)}

(the division by :math:`T_{tot}` is only applied where the numerator is positive).

Aerosol retrieval
-----------------

:math:`\tau_{ref}` is the largest AOT that does not make the reflectance of too many water pixels
negative. For a set of tested values, the number of pixels with a negative minimum reflectance is

.. math::

   N(\tau_{ref}) = \#\left\{ \mathbf{x} \in V \;:\; \min_\lambda R(\lambda, \mathbf{x}; \tau_{ref}) < 0 \right\},
   \qquad V = \left\{ \mathbf{x} : \min_\lambda R(\lambda, \mathbf{x}; 0) > t_0 \right\}

and the retrieved value is

.. math::

   \hat\tau_{ref} = \max\left\{ \tau_{ref} \;:\; N(\tau_{ref}) < \eta\, \max\left(N_0, \max_\tau N(\tau)\right) \right\}

with :math:`\eta` = ``neg_pix_rel_thresh``. It is first bounded with the SWIR bands
(`Kernel.aerosol_swir_chunk`), then retrieved from the visible-NIR bands, one value per block of the
coarse raster (`Kernel.aerosol_visible_chunk`).

Processing sequence
===================

Typical sequence (``prod`` = product with ``prod.raster``, ``cams`` = CAMS ancillary data)::

    kernel = Kernel(prod, aero_lut, trans_lut, cams)
    kernel.set_gas_transmittance()
    kernel.lut_preparation()
    coarse = kernel.get_coarse_masked_raster(mask=mask)         # water pixels, coarser grid

    wind = kernel.rugosity_est_chunk(coarse)                     # wind speed from sunglint
    kernel.wind_img = kernel.smoothing(wind)                     # (small) wind map, used by next steps
    kernel.aerosol_swir_chunk(coarse)                            # sets kernel.aot_ref_max
    aot = kernel.aerosol_visible_chunk(coarse_block)             # one value per block, see below
    kernel.final_process(prod.raster, aot_ref_raster)            # lazy: kernel.xres (Rrs, BRDFg, aot550)

Dask design
===========

* **Block functions are pure.** ``*_chunk`` methods and ``final_process_chunk`` only *read* the
  kernel (LUTs, ancillary data, settings) and return their result: nothing is written to ``self`` while a
  block is processed (the previous version stored ``Tdir``, ``_p_slope_``, ``aot_ref_max``... on the
  instance, which is a race between threads and keeps the arrays of the last block alive). They can be
  used in ``xr.map_blocks`` / ``dask.delayed`` with the threaded or the distributed scheduler.
* **Blocks are loaded once.** A block may be lazy (dask-backed); each function loads only the variables and
  bands it needs a single time, instead of re-evaluating the graph at every ``.values`` of an
  ``aot_ref`` loop.
* **Small data are loaded once.** Look-up tables, coarse gaseous transmittance, surface pressure and CAMS
  AOD are small: they are computed once (``_load_small``) so that blocks only interpolate them. Their
  size limit is ``_SMALL_ARRAY_BYTES``.
* **Nothing full-size is allocated.** ``final_process`` is lazy: it returns a dask-backed Dataset built
  with ``xr.map_blocks``, evaluated chunk by chunk when it is written (e.g. by ``L2aProduct``), so the
  full ``Rrs`` raster is never in memory. The ``aot_ref`` loops keep a running pixel count instead of one
  full raster per ``aot_ref``.
* **Numba:** the kernels of ``Rasterization`` are called from several dask threads. If they are compiled
  with ``parallel=True``, the default ``workqueue`` threading layer aborts on concurrent calls: use
  ``NUMBA_THREADING_LAYER=omp`` (or ``tbb``), or a process-based scheduler.

Chunking: ``final_process(..., chunk=512)`` divides the 1024 x 1024 shards of ``L2aProduct.export_to_zarr``,
which merges four blocks per shard at export time without recomputing them.
'''
import copy
import logging
import types
import warnings

import numpy as np
import xarray as xr

# keep attributes through operation on xarray objects
xr.set_options(keep_attrs=True)
import rioxarray  # noqa: F401  (registers the .rio accessor used by the callers)
import yaml

from importlib.resources import files
from numba import njit
from scipy import ndimage

from .auxdata import AuxData
from .utils import Rasterization, GaseousTransmittance, Misc

# --------------------------------------------------
# get path of other files as indicated in config.yml
# --------------------------------------------------
configfile = files(__package__) / 'config.yml'
with open(configfile, 'r') as file:
    config = yaml.safe_load(file)

AEROSOL_MODELS = config['settings']['aerosol_models']
AEROSOL_COMBINATION = config['settings']['aerosol_combination']

# arrays smaller than this (bytes) are loaded in memory once instead of being recomputed by every block
_SMALL_ARRAY_BYTES = 512 * 2 ** 20


def _load_small(obj, max_bytes=_SMALL_ARRAY_BYTES):
    '''
    Load a small dask-backed xarray object in memory (once); other objects are returned unchanged.

    :param obj: DataArray / Dataset (anything else is returned as is)
    :param max_bytes: objects larger than this stay lazy (they are not "small" data)
    :return: numpy-backed object if `obj` was a small dask-backed one, else `obj`
    '''
    if isinstance(obj, (xr.DataArray, xr.Dataset)) and obj.chunks and obj.nbytes <= max_bytes:
        return obj.compute()
    return obj


def _unique_values(*arrays):
    '''
    Sorted unique values of several DataArrays in one pass over their data (dask-aware: only the unique
    values of each chunk are gathered, the arrays are never fully in memory).

    :param arrays: DataArrays (dask or numpy backed)
    :return: tuple of numpy arrays, one per input
    '''
    try:
        import dask
        import dask.array as dsa
    except ImportError:
        return tuple(np.unique(np.asarray(a)) for a in arrays)
    lazy = [dsa.unique(a.data) if isinstance(a.data, dsa.Array) else np.unique(a.data) for a in arrays]
    return dask.compute(*lazy)


class _NegativeCounter:
    '''
    Running count of pixels with negative reflectance for increasing aerosol optical thickness.

    The AOT retrievals look for the largest ``aot_ref`` for which only few pixels get a negative minimum
    reflectance. Instead of stacking one full raster per ``aot_ref``, only the count is kept:

    * the first ``add`` (``aot_ref = 0``) defines the valid pixels: ``metric > base_threshold``
      (this removes pixels that are already negative after the Rayleigh correction, e.g. shadows);
    * each ``add`` counts the valid pixels whose metric is negative.
    '''

    def __init__(self, base_threshold):
        '''
        :param base_threshold: a pixel is valid if its metric at ``aot_ref = 0`` is above this value
        '''
        self.base_threshold = base_threshold
        self._valid = None
        self._counts = []

    def add(self, reflectance, nwl=None):
        '''
        Add the reflectance of the next ``aot_ref``.

        :param reflectance: array (wl, y, x)
        :param nwl: only the first `nwl` bands are used (None: all)
        '''
        values = np.asarray(reflectance)
        if nwl is not None:
            values = values[:nwl]
        with warnings.catch_warnings(), np.errstate(invalid='ignore'):
            warnings.simplefilter('ignore', RuntimeWarning)
            metric = np.nanmin(values, axis=0)  # minimum over wavelengths, NaN skipped
            if self._valid is None:
                self._valid = metric > self.base_threshold
            self._counts.append(int(np.count_nonzero(self._valid & (metric < 0))))

    @property
    def counts(self):
        '''Number of negative pixels per added ``aot_ref`` (int64 array).'''
        return np.asarray(self._counts, dtype=np.int64)


def _aot_ref_candidates(aot_refs, counts, previous_max, rel_thresh):
    '''
    ``aot_ref`` values with few negative pixels, in the order of `aot_refs`.

    :param aot_refs: tested ``aot_ref`` values
    :param counts: number of negative pixels for each of them (see `_NegativeCounter`)
    :param previous_max: lower bound of the reference number of negative pixels
    :param rel_thresh: accepted fraction of the reference number (max of `previous_max` and `counts`)
    :return: array of the ``aot_ref`` whose count is below ``rel_thresh * reference``
    '''
    neg_pix_max = np.max([previous_max, counts.max()])
    return np.asarray(aot_refs)[counts < neg_pix_max * rel_thresh]


@njit
def _find_min_cost(cost, values):
    '''
    Value with the smallest cost along the first axis, pixel by pixel.

    Defined once at module level: a function defined inside a method is compiled again by numba at every call.

    :param cost: array (n_values, ny, nx)
    :param values: candidate values, length n_values
    :return: ``(best_value, best_cost)``, two arrays (ny, nx); NaN where the cost is NaN for all values
    '''
    n, ny, nx = cost.shape
    best_value = np.full((ny, nx), np.nan)
    best_cost = np.full((ny, nx), np.nan)
    for iy in range(ny):
        for ix in range(nx):
            if np.isnan(cost[:, iy, ix]).all():
                continue
            ind = np.argmin(cost[:, iy, ix])
            best_value[iy, ix] = values[ind]
            best_cost[iy, ix] = cost[ind, iy, ix]
    return best_value, best_cost


class Kernel():
    '''
    Retrieval of wind speed and aerosol optical thickness over water, and computation of Rrs.

    Parameters
    ----------
    prod : object
        Product with ``prod.raster`` (Dataset with ``bands``, ``sza``, ``vza``, ``raa``, ``wl_true``),
        ``prod.sensor`` (string, ``'S2'`` for Sentinel-2) and ``prod._type`` (dtype of the output arrays).
    aero_lut : xarray.Dataset
        Aerosol look-up table (``I``, ``aot``; dimensions ``wind``, ``model``, ``aot_ref``, ``sza``, ``vza``,
        ``azi``, ``wl``).
    trans_lut : xarray.Dataset
        Transmittance look-up table (``Ttot_Ed``...).
    cams : object
        CAMS ancillary data: ``cams.raster.sp`` (surface pressure) and ``cams.cams_aod``.

    Notes
    -----
    * The statistics of the CAMS AOD at 550 nm (``aot_ref_cams_max`` / ``_mean`` / ``_min``) are plain floats.
    * Small arrays (surface pressure, CAMS AOD, LUTs, gaseous transmittance) are loaded in memory once and
      cached. Replace ``cams`` or the LUTs by creating a new ``Kernel``.
    * Block functions (``*_chunk``) do not modify the instance, except ``aerosol_swir_chunk`` which sets
      the global ``aot_ref_max`` (it works on the whole coarse raster).
    '''

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
        self.aot_ref_cams_max = float(aot_ref_cams.max())
        self.aot_ref_cams_mean = float(aot_ref_cams.mean())
        self.aot_ref_cams_min = float(aot_ref_cams.min())

        self.wl_true = prod.raster.wl_true
        self.pressure_ref = 101320.

        self.neg_pix_max = 0

        if 'S2' in prod.sensor:
            self.monoview = False
        else:
            self.monoview = True
        self._R_ = Rasterization(monoview=self.monoview)
        self._cache = {}

    # ------------------------------------------------------------------
    # helpers shared by the block functions
    # ------------------------------------------------------------------
    def _cached(self, key, factory):
        '''Value computed once by `factory` and kept (concurrent first calls just compute it twice).'''
        try:
            return self._cache[key]
        except KeyError:
            value = self._cache[key] = factory()
            return value

    @property
    def _surface_pressure(self):
        '''CAMS surface pressure, in memory (coarse grid).'''
        return self._cached('sp', lambda: _load_small(self.cams.raster.sp))

    @property
    def _cams_aod550(self):
        '''CAMS AOD at 550 nm, in memory (coarse grid).'''
        return self._cached('aod550', lambda: _load_small(self.cams.cams_aod.interp(wl=550, method='quadratic')))

    @staticmethod
    def _load_block(raster, variables):
        '''
        Load the needed variables of a block in memory, once.

        :param raster: Dataset (dask or numpy backed)
        :param variables: names of the variables to keep (the ones missing in `raster` are ignored)
        :return: numpy-backed Dataset restricted to those variables
        '''
        names = [v for v in variables if v in raster.data_vars]
        return raster[names].compute()

    # ------------------------------------------------------------------
    # preparation (whole image / look-up tables)
    # ------------------------------------------------------------------
    def coarsening_factor(self, resolution):
        r'''
        Number of pixels per coarse pixel, along x and y, to reach the target resolution.

        With :math:`\Delta x`, :math:`\Delta y` the pixel size of the raster and :math:`r` the target
        resolution:

        .. math::

           f_x = \frac{r}{\Delta x}, \qquad f_y = \frac{r}{\Delta y}, \qquad f_x, f_y \in \mathbb{N}^*

        e.g. a 20 m image gives :math:`f = 3, 6, 12` for :math:`r` = 60, 120, 240 m.

        :param resolution: target pixel size, in the unit of the ``x`` / ``y`` coordinates (meter),
            e.g. 60, 120, 240...; must be a multiple of the pixel size of the raster
        :return: dict ``{'x': factor, 'y': factor}`` for ``xarray.Dataset.coarsen``
        '''
        factors = {}
        for dim in ('x', 'y'):
            pixel_size = abs(float(self.raster[dim][1] - self.raster[dim][0]))
            factor = resolution / pixel_size
            if round(factor) < 1 or not np.isclose(factor, round(factor)):
                raise ValueError(f'coarse resolution {resolution} is not a multiple of the pixel size '
                                 f'{pixel_size} along {dim}')
            factors[dim] = int(round(factor))
        return factors

    def get_coarse_masked_raster(self,
                                 variables=['sza', 'vza', 'raa', 'bands'],

                                 mask=None,
                                 resolution=240,
                                 min_pixel_prop=0.5):
        r'''
        Coarse-resolution raster of the masked image (lazy if the raster is dask-backed).

        The pixels with ``mask != 0`` are discarded, the others are averaged by blocks of coarse pixels of
        size ``resolution`` (incomplete blocks at the edges are padded). Coarse pixels made of too few
        valid pixels are set to NaN.

        For the coarse pixel :math:`B` made of the :math:`f_x \times f_y` pixels of the block
        (see `coarsening_factor`) and :math:`W_B \subset B` its valid (water) pixels:

        .. math::

           \bar X_B = \frac{1}{N_B} \sum_{\mathbf{x} \in W_B} X(\mathbf{x}), \qquad N_B = \# W_B

        for each variable :math:`X`, and :math:`B` is kept only if

        .. math::

           N_B > p \, \max_{B'} N_{B'}

        with :math:`p` = ``min_pixel_prop``. :math:`N_B` is returned as ``water_pixel_number``.

        :param variables: variables of the raster to coarsen
        :param mask: raster of the same grid, 0 for the pixels to keep (None: keep all)
        :param resolution: pixel size (meter) of the coarse raster, e.g. 60, 120, 240...; must be a
            multiple of the pixel size of the raster (see ``coarsening_factor``)
        :param min_pixel_prop: a coarse pixel is kept if its number of valid pixels is above this proportion
            of the largest number of valid pixels of a coarse pixel
        :return: coarse Dataset with the variables and ``water_pixel_number`` (also kept in
            ``self.coarse_masked_raster``, unfiltered); the resolution is kept in ``self.coarse_resolution``
            and in the ``resolution`` attribute
        '''
        window = self.coarsening_factor(resolution)
        self.coarse_resolution = resolution
        if mask is not None:
            raster = self.raster[variables].where(mask == 0)
        else:
            raster = self.raster[variables]
        self.coarse_masked_raster = raster.coarsen(window, boundary="pad").mean()
        self.coarse_masked_raster['water_pixel_number'] = raster['sza']. \
            coarsen(window, boundary="pad").count()
        self.coarse_masked_raster.attrs['resolution'] = resolution
        min_pixel_number = min_pixel_prop * np.max(self.coarse_masked_raster['water_pixel_number'])
        return self.coarse_masked_raster.where(self.coarse_masked_raster.water_pixel_number > min_pixel_number)

    def lut_preparation(self,
                        wind=2,
                        ang_resol={'sza': 1, 'vza': 1, 'raa': 0},
                        aot_refs=np.linspace(0.0, 0.8, 25),
                        aerosol_combination=AEROSOL_COMBINATION):
        r'''
        Reduce the look-up tables to the geometries of the image and to the aerosol mixture.

        The aerosol models :math:`m` are mixed with the weights :math:`w_m` (`aerosol_combination`) and the
        LUT intensities are converted into reflectance:

        .. math::

           R_{diff}(\lambda; \tau_{ref}) = \frac{1}{\cos\theta_s} \sum_m w_m\, I_m(\lambda; \tau_{ref}),
           \qquad
           R_{ray}(\lambda) = \frac{I_0(\lambda; \tau_{ref} = 0)}{\cos\theta_s}

           \tau_a(\lambda; \tau_{ref}) = \sum_m w_m\, \tau_{a,m}(\lambda; \tau_{ref}),
           \qquad
           T_{\downarrow}(\lambda; \theta, \tau_{ref}) = \sum_m w_m\, T_{\downarrow,m}(\lambda; \theta, \tau_{ref})

        (:math:`I_0`: first model of the LUT, used for the Rayleigh-only reflectance). All the tables are
        then interpolated on `aot_refs` (quadratic) and on the band wavelengths :math:`\lambda`.

        The angles are rounded (`ang_resol`) and only their unique values are kept in the LUTs, which
        speeds up the later interpolations. The unique values are found in one pass over the angle rasters
        without loading them (dask). The resulting tables are loaded in memory once, so the blocks do not
        recompute them.

        Sets ``trans_aero_lut``, ``Rray``, ``Rdiff_lut``, ``aot_lut``, ``szas``, ``vzas``, ``azis``,
        ``aot_refs``, ``sunglint_eps``, ``rot``, and adds ``sza_trunc``, ``vza_trunc``, ``raa_trunc`` to
        ``prod.raster``.

        :param wind: wind speed (m/s) of the LUT
        :param ang_resol: number of decimals kept for each angle
        :param aot_refs: AOT at the reference wavelength of the interpolated LUTs
        :param aerosol_combination: weights of the aerosol models (array, list or DataArray along ``model``)
        '''

        logging.info('LUT preparation')

        if isinstance(aerosol_combination, (list, np.ndarray)):
            aerosol_combination = xr.DataArray(aerosol_combination, coords={'model': AEROSOL_MODELS})
        elif not isinstance(aerosol_combination, xr.DataArray):
            logging.error("error in aerosol_combination parameter, should be numpy or xr.DataArray")
            return

        self.aerosol_combination = aerosol_combination

        aero_lut = self.aero_lut.sel(wind=wind, method='nearest')
        trans_lut = self.trans_lut.sel(wind=wind, method='nearest')

        # -----------------------------------------------------------------
        # Simplify (truncation) of angles values
        # to speed up further interpolation based on nearest neighbors
        # -----------------------------------------------------------------
        for param in ['sza', 'vza', 'raa']:
            self.prod.raster[param + '_trunc'] = self.prod.raster[param].round(ang_resol[param])

        # unique values of the three angle rasters, one pass over the data; cached because they only depend
        # on the angles and their rounding, so the second call (refined LUT) does not read the rasters again
        sza_, vza_, azi_ = self._cached(('angles', tuple(sorted(ang_resol.items()))),
                                        lambda: _unique_values(self.prod.raster.sza_trunc,
                                                               self.prod.raster.vza_trunc,
                                                               (180. - self.prod.raster.raa_trunc) % 360))

        sza_ = sza_[~np.isnan(sza_)]
        vza_ = vza_[~np.isnan(vza_)]
        azi_ = azi_[~np.isnan(azi_)]

        # ----------------------------------------
        # interpolation transmittance
        # ----------------------------------------
        self.trans_aero_lut = trans_lut.interp(sza=[*sza_, *vza_])
        self.trans_aero_lut = (self.trans_aero_lut * aerosol_combination).sum('model')
        self.trans_aero_lut = self.trans_aero_lut.interp(aot_ref=aot_refs, method='quadratic')
        self.trans_aero_lut = self.trans_aero_lut.interp(wl=self.wl_true, method='quadratic')
        # clean up the xarray DataArray object:
        self.trans_aero_lut = self.trans_aero_lut.squeeze().reset_coords(drop=True)

        # ----------------------------------------
        # interpolation Rayleigh
        # ----------------------------------------
        self.Rray = aero_lut.I.isel(model=0).sel(aot_ref=0, method='nearest')
        self.Rray = self.Rray.interp(sza=sza_, vza=vza_).interp(azi=azi_)
        self.Rray = (self.Rray.interp(wl=self.wl_true, method='quadratic') /
                     np.cos(np.radians(self.Rray.sza)))

        # ----------------------------------------
        # interpolation atmo diffuse light
        # ----------------------------------------
        self.Rdiff_lut = aero_lut.I.interp(sza=sza_, vza=vza_).interp(azi=azi_)
        self.Rdiff_lut = (self.Rdiff_lut * aerosol_combination).sum('model')
        self.Rdiff_lut = self.Rdiff_lut.interp(aot_ref=aot_refs, method='quadratic')
        self.Rdiff_lut = (self.Rdiff_lut.interp(wl=self.wl_true, method='quadratic') /
                          np.cos(np.radians(self.Rdiff_lut.sza)))

        # ----------------------------------------
        # interpolation aerosol optical thickness
        # ----------------------------------------
        self.aot_lut = (aero_lut.aot * aerosol_combination).sum('model')
        self.aot_lut = self.aot_lut.interp(aot_ref=aot_refs,
                                           method='quadratic'
                                           ).interp(wl=self.wl_true, method='cubic')

        # the tables are used by every block: compute them once if they are lazy
        self.trans_aero_lut = _load_small(self.trans_aero_lut)
        self.Rray = _load_small(self.Rray)
        self.Rdiff_lut = _load_small(self.Rdiff_lut)
        self.aot_lut = _load_small(self.aot_lut)

        self.szas = self.Rdiff_lut.sza.values
        self.vzas = self.Rdiff_lut.vza.values
        self.azis = self.Rdiff_lut.azi.values
        self.aot_refs = self.Rdiff_lut.aot_ref.values

        _auxdata = AuxData(wl=self.wl_true)  # wl=masked.wl)
        self.sunglint_eps = _load_small(_auxdata.sunglint_eps)  # ['mean'].interp(wl=wl_true)
        self.rot = _load_small(_auxdata.rot)

    def mu_N(self,
             sza,
             vza,
             azi,
             monoview=False):
        r'''
        Compute the normal angle to wave slopes that produce sunglint.
        Warning: azi: azimuth in rad for convention azi=180 when sun-sensenor in oppositioon

        The facet reflecting the sun towards the sensor bisects the incident and reflected directions;
        with the scattering angle :math:`\alpha` between them:

        .. math::

           \cos\alpha = \cos\phi \sin\theta_v \sin\theta_s + \cos\theta_v \cos\theta_s,
           \qquad
           \cos\theta_N = \frac{\cos\theta_s + \cos\theta_v}{\sqrt{2\,(1 + \cos\alpha)}}

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
        r'''
        Compute propability of wave slopes producing sunglint.

        Isotropic Cox-Munk (1954) distribution of the slope of the facet normal :math:`\theta_N`
        (see `mu_N`), divided by :math:`\cos^4\theta_N`:

        .. math::

           P_{slope} = \frac{1}{\pi \sigma^2}\,
           \exp\left(-\frac{\tan^2\theta_N}{\sigma^2}\right) \frac{1}{\cos^4\theta_N},
           \qquad \sigma^2 = \frac{U + 0.586}{195.3}

        Warning: the correction is currently disabled (see the TODO below): the function returns 1
        everywhere, and NaN where the angles or `sigma2` are NaN (this NaN pattern is used downstream).
        The slope probability is still evaluated, so while the correction stays disabled the result does not
        depend on `sigma2` for finite values (it could be computed once outside the wind loop of
        ``rugosity_est_chunk``).

        :param sza: solar zenith angle in degree
        :param vza: viewing zenith angle in degree
        :param azi: relative azimuth bewteen sun and sensor
        :param sigma2: mean square slope of the wave slope distribution (scalar or array y, x)
        :param monoview: Set sensor viewing  configuration:
            - monoview = True : same viewing angles for all the spectral bands
            - monoview = False : viewing angles depend on spectral band (e.g. Sentinel-2 images)

        :return: DataArray of ones (y, x) if monoview, else (wl, y, x)
        '''

        cosN = self.mu_N(sza, vza, azi, monoview=monoview)
        thetaN = np.arccos(cosN)
        # stats == 'cm_iso':
        # TODO check consistency between sigma2 and formulation
        # Pdist_ = 1. / (np.pi *2.* sigma2) * np.exp(-1./2 * np.tan(thetaN) ** 2 / sigma2)
        xp_slope = 1. / (np.pi * sigma2) * np.exp(- np.tan(thetaN) ** 2 / sigma2) / cosN ** 4
        if monoview:
            xp_slope = xp_slope.transpose('y', 'x')
        else:
            xp_slope = xp_slope.transpose('wl', 'y', 'x')
            # TODO Warning revise equation issue with azimuth: correction discarded!!
        return xp_slope * 0 + 1

    def set_gas_transmittance(self):
        r'''
        Gaseous transmittance (H2O, CH4) of the low altitude gases, on the coarse grid.

        Sets ``Tg_diff_raster`` (diffuse light, absorption-scattering coefficients 0.5 / 0.65) and
        ``Tg_raster`` (coefficients 1), dimensions (wl, y, x). They are loaded in memory once (coarse data):
        the blocks interpolate them to their own grid.

        With :math:`T_{gas}(\lambda; c)` the band transmittance of a gas for the column :math:`c` and
        :math:`C_{gas}` the CAMS total column:

        .. math::

           T_g(\lambda) = \prod_{gas \in \{H_2O,\, CH_4\}} T_{gas}(\lambda;\, C_{gas}),
           \qquad
           T_g^{diff}(\lambda) = T_{H_2O}(\lambda;\, 0.5\, C_{H_2O})\; T_{CH_4}(\lambda;\, 0.65\, C_{CH_4})

        The coefficients :math:`k_{gas} < 1` of :math:`T_g^{diff}` account for the gases located below the
        scattering layer, which absorb only part of the diffuse light path (6S approach).
        '''
        logging.info('set gaseous transmittance')
        # -------------------------------------------------------------
        # SET GASEOUS TRANSMITTANCE FOR LOW ALTITUDE GASES
        # -------------------------------------------------------------
        gases = ['h2o', 'ch4']
        gas_trans = GaseousTransmittance(self.prod, self.cams)
        # set total transmittance values
        gas_trans.coef_abs_scat['h2o'] = 0.5
        gas_trans.coef_abs_scat['ch4'] = 0.65
        self.Tg_diff_raster = _load_small(
            gas_trans.get_gaseous_transmittance(gases=gases, background=False).transpose("wl", "y", "x"))

        # set total transmittance values
        for gas in gases:
            gas_trans.coef_abs_scat[gas] = 1
        self.Tg_raster = _load_small(
            gas_trans.get_gaseous_transmittance(gases=gases, background=False).transpose("wl", "y", "x"))

    # ------------------------------------------------------------------
    # block functions
    # ------------------------------------------------------------------
    def rugosity_est_chunk(self,
                           raster,
                           aot_ref=0.1,
                           winds=np.linspace(0.2, 15, 25),
                           idx_nir=[6, 7, 8, 10],
                           chunk=100000
                           ):
        r'''
        Wind speed retrieval from the sunglint, on one block (use with ``xr.map_blocks``).

        For each candidate wind speed, the glint-free reflectance is extrapolated between the NIR bands
        (parallax / isotropic Cox-Munk model, see `p_slope`); the retrieved wind is the one that minimises
        the squared residual. Needs ``set_gas_transmittance`` and ``lut_preparation``.

        With the bands :math:`\lambda_0, \lambda_1, \lambda_2` = ``idx_nir[0:3]`` and :math:`R_{corr}`
        computed for ``aot_ref`` (see the module documentation), the glint is normalised to the geometry
        of :math:`\lambda_2`:

        .. math::

           R_{geom}(\lambda_i; U) = \frac{P_{slope}(\lambda_2; U)}{T_g(\lambda_i)\, P_{slope}(\lambda_i; U)}\,
           R_{corr}(\lambda_i)

        and the residual is the difference between :math:`R_{geom}(\lambda_1)` and its linear
        interpolation between :math:`\lambda_0` and :math:`\lambda_2`:

        .. math::

           \Delta R(U) = R_{geom}(\lambda_0) + \left(R_{geom}(\lambda_2) - R_{geom}(\lambda_0)\right)
           \frac{\lambda_1 - \lambda_0}{\lambda_2 - \lambda_0} - R_{geom}(\lambda_1),
           \qquad
           \hat U = \underset{U}{\operatorname{argmin}}\ \Delta R(U)^2

        Warning: :math:`U` only acts through :math:`P_{slope}`, which is currently disabled (see
        `p_slope`): :math:`\Delta R` is then the same for all the candidates and the first wind of `winds`
        is returned.

        :param raster: Dataset block (lazy or not) with ``bands``, ``sza``, ``raa``, ``vza`` and the
            coordinates ``wl``, ``x``, ``y``
        :param aot_ref: AOT at the reference wavelength used for the diffuse light correction
        :param winds: candidate wind speeds (m/s)
        :param idx_nir: indices of the (at least four) NIR / SWIR bands used
        :param chunk: unused, kept for compatibility
        :return: Dataset (y, x) with ``wind`` (m/s, float64), ``cost`` (float64) and ``brdf_g`` (the SWIR
            reflectance after the diffuse light correction), or None for an empty block
        '''

        monoview = self.monoview
        _R_ = self._R_

        raster = self._load_block(raster.isel(wl=idx_nir), ('bands', 'sza', 'raa', 'vza'))
        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        if Ny == 0 or Nx == 0:
            return

        _wl_true = self.wl_true[idx_nir]
        Rdiff_lut_ = self.Rdiff_lut.isel(wl=idx_nir)

        _pressure = self._surface_pressure.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.isel(wl=idx_nir).interp(x=raster.x, y=raster.y)
        _Tg_raster = self.Tg_raster.isel(wl=idx_nir).interp(x=raster.x, y=raster.y)

        # subsetting
        _sza = raster.sza
        _raa = raster.raa
        _vza = raster.vza
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

        best_wind, best_cost = _find_min_cost(delta_R.values ** 2, np.asarray(winds, dtype=np.float64))

        return xr.Dataset(dict(wind=(["y", "x"], best_wind),
                               cost=(["y", "x"], best_cost),
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
                           chunk=100000,
                           keep_xres=False
                           ):
        r'''
        Maximum AOT compatible with the SWIR bands, on the whole (coarse masked) raster.

        For each tested ``aot_ref`` the SWIR reflectance is corrected for the diffuse light of the
        atmosphere; ``aot_ref_max`` is the largest tested value for which fewer than
        ``neg_pix_rel_thresh`` of the maximum number of pixels get a negative reflectance (among the pixels
        positive at ``aot_ref = 0``). The counts are accumulated on the fly: only one reflectance array is
        in memory at a time.

        With the two SWIR bands :math:`\Lambda_s` = ``idx_swir``:

        .. math::

           R_{corr}(\lambda; \tau_{ref}) = R_{TOA}(\lambda) - T_g^{diff}(\lambda)\, \frac{P}{P_0}\,
           R_{diff}(\lambda; \tau_{ref}), \qquad \lambda \in \Lambda_s

           N(\tau_{ref}) = \#\left\{\mathbf{x} \in V : \min_{\lambda \in \Lambda_s} R_{corr}(\lambda, \mathbf{x};
           \tau_{ref}) < 0\right\},
           \qquad V = \left\{\mathbf{x} : \min_{\lambda \in \Lambda_s} R_{corr}(\lambda, \mathbf{x}; 0) > 0\right\}

           \tau_{max} = \max\left\{\tau_{ref} : N(\tau_{ref}) < \eta\, \max_\tau N(\tau)\right\}

        with :math:`\eta` = ``neg_pix_rel_thresh``. The water-leaving signal is negligible in the SWIR: an
        AOT above :math:`\tau_{max}` would remove more than the atmospheric signal.

        Sets ``self.aot_ref_max``; if the result is outside ``]0, aot_ref_lim[`` it falls back on
        ``min(aot_ref_cams_max, aot_ref_lim)``; with no candidate it starts from ``0.8 * aot_refs[-1]``.

        :param raster: Dataset (coarse masked raster, lazy or not) with ``bands``, ``sza``, ``raa``, ``vza``
        :param aot_refs: tested AOT values at the reference wavelength, increasing, starting at 0
        :param idx_swir: indices of the two SWIR bands
        :param neg_pix_rel_thresh: accepted fraction of negative pixels
        :param chunk: unused, kept for compatibility
        :param keep_xres: keep the reflectance of every ``aot_ref`` in ``self.xres`` (debugging; one full
            raster per ``aot_ref``)
        :return: ``aot_ref_max``, or None for an empty raster
        '''

        _R_ = self._R_

        raster = self._load_block(raster.isel(wl=idx_swir), ('bands', 'sza', 'raa', 'vza'))
        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        if Ny == 0 or Nx == 0:
            return

        Rdiff_lut_ = self.Rdiff_lut.isel(wl=idx_swir)

        _pressure = self._surface_pressure.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.isel(wl=idx_swir).interp(x=raster.x, y=raster.y)
        _Tg_abs_diff = _Tg_diff_raster.values

        # subsetting (loaded once, the aot_ref loop only reuses them)
        _sza = raster.sza.values
        _raa = raster.raa
        _vza = raster.vza.values
        _azi = ((180. - _raa) % 360).values

        _pressure_ = _pressure / self.pressure_ref
        band = _band_rad.values

        counter = _NegativeCounter(base_threshold=0.)
        xres = []
        for _aot_ref in aot_refs:
            # get LUT values
            _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, _sza,
                                              self.vzas, _vza,
                                              self.azis, _azi,
                                              Nwl, Ny, Nx, Rdiff_lut_.sel(aot_ref=_aot_ref, method='nearest').values)

            _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

            #  correction for diffuse light
            Rcorr = band - _Rdiff
            counter.add(Rcorr)
            if keep_xres:
                xres.append(xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rcorr)),
                                       coords=dict(wl=raster.wl, x=raster.x, y=raster.y, aot_ref=_aot_ref)))

        if keep_xres:
            self.xres = xr.concat(xres, dim='aot_ref')

        # get aot550 max from which negative Rrs appear
        aot_ref_maxs = _aot_ref_candidates(aot_refs, counter.counts, 0, neg_pix_rel_thresh)
        if len(aot_ref_maxs) > 0:
            aot_ref_max = aot_ref_maxs[-1]
        else:
            aot_ref_max = aot_refs[-1] * 0.8  # self.aot_ref_cams_max

        # to stay within LUT aot range
        if not ((aot_ref_max > 0) and (aot_ref_max < self.aot_ref_lim)):
            aot_ref_max = np.min([self.aot_ref_cams_max, self.aot_ref_lim])

        self.aot_ref_max = aot_ref_max
        return aot_ref_max

    def aerosol_visible_chunk(self,
                              raster,
                              # aot_refs=[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8],
                              iwl_swir=[-2, -1],
                              neg_pix_rel_thresh=0.05,
                              chunk=100000,
                              keep_xres=False,
                              update_state=False
                              ):
        r'''
        AOT retrieval from the visible-NIR bands, on one block (use with ``xr.map_blocks`` / ``dask.delayed``).

        For 21 values of ``aot_ref`` between 0 and ``min(1.02 * aot_ref_max, aot_ref_lim)``, the reflectance
        is corrected for the diffuse light and the sunglint, and the pixels whose minimum reflectance over
        the first 7 bands is negative are counted (among the pixels above -0.0003 at ``aot_ref = 0``, which
        removes e.g. shadows). The retrieved AOT is the largest value with fewer negative pixels than
        ``neg_pix_rel_thresh`` of the maximum count; without such value it is the mean CAMS AOD over the
        block. The counts are accumulated on the fly (one reflectance array in memory at a time).

        For each tested :math:`\tau_{ref}`, with :math:`R_{corr}` and :math:`T_{dir}` as in the module
        documentation (here :math:`P/P_0` is averaged over the block in :math:`T_{dir}`), the sunglint is
        estimated in both SWIR bands :math:`\lambda_s \in \Lambda_s` = ``iwl_swir`` and normalised to the
        last band :math:`\lambda_r`:

        .. math::

           R_f(\lambda_s) = \frac{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)\, R_{corr}(\lambda_s)}
           {T_g(\lambda_s)\, T_{dir}(\lambda_s)\, \varepsilon(\lambda_s)\, P_{slope}(\lambda_s)},
           \qquad
           R_f = \min_{\lambda_s \in \Lambda_s} \left|R_f(\lambda_s)\right|

           R_g(\lambda) = T_g(\lambda)\, T_{dir}(\lambda)\, \varepsilon(\lambda)\, P_{slope}(\lambda)\,
           \frac{R_f}{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)},
           \qquad
           R_{rs}^{TOA}(\lambda; \tau_{ref}) = \frac{R_{corr}(\lambda) - R_g(\lambda)}{\pi}

        (for a mono-view sensor, :math:`R_f(\lambda_s)` is not multiplied by
        :math:`\varepsilon(\lambda_r) P_{slope}(\lambda_r)`). The negative pixels are counted on the first
        seven bands :math:`\Lambda_7`:

        .. math::

           N(\tau_{ref}) = \#\left\{\mathbf{x} \in V : \min_{\lambda \in \Lambda_7}
           R_{rs}^{TOA}(\lambda, \mathbf{x}; \tau_{ref}) < 0\right\},
           \qquad V = \left\{\mathbf{x} : \min_{\lambda \in \Lambda_7} R_{rs}^{TOA}(\lambda, \mathbf{x}; 0)
           > -0.0003\right\}

           \hat\tau_{ref} = \max\left\{\tau_{ref} : N(\tau_{ref}) < \eta\,
           \max\left(N_0, \max_\tau N(\tau)\right)\right\}

        with :math:`\eta` = ``neg_pix_rel_thresh`` and :math:`N_0` = ``self.neg_pix_max``.

        The block only reads ``self.aot_ref_max`` (set by ``aerosol_swir_chunk``) and ``self.neg_pix_max``;
        it does not modify the instance, so blocks are independent of each other and of the order or the
        threads they run in (the previous version overwrote ``aot_ref_max`` at the end of every block).

        :param raster: Dataset block (lazy or not) with ``bands``, ``sza``, ``raa``, ``vza``,
            ``water_pixel_number`` and the coordinates ``wl``, ``x``, ``y``
        :param iwl_swir: indices of the two SWIR bands used for the sunglint
        :param neg_pix_rel_thresh: accepted fraction of negative pixels
        :param chunk: unused, kept for compatibility
        :param keep_xres: keep the reflectance of every ``aot_ref`` in ``self.xres`` (debugging; 21 full
            blocks, not thread safe)
        :param update_state: store the retrieved AOT in ``self.aot_ref_max`` (only meaningful when the
            function is called once on the whole raster; not thread safe)
        :return: Dataset with one pixel (y, x = block centre): ``aot_ref`` and ``number_of_pixel`` (number
            of water pixels of the block), or None for an empty block
        '''
        ##############################################
        # Estimate AOT from visible-NIR
        ##############################################
        monoview = self.monoview
        _R_ = self._R_

        raster = self._load_block(raster, ('bands', 'sza', 'raa', 'vza', 'water_pixel_number'))
        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        if Ny == 0 or Nx == 0:
            return

        # settings read once: they must not change while the block is processed
        aot_ref_max_in = self.aot_ref_max
        neg_pix_max_in = self.neg_pix_max
        sunglint_eps = self.sunglint_eps.values
        rot = self.rot.values

        _pressure = self._surface_pressure.interp(x=raster.x, y=raster.y).values
        _Tg_diff_raster = self.Tg_diff_raster.interp(x=raster.x, y=raster.y)
        _Tg_abs_diff = _Tg_diff_raster.values
        _Tg_raster = self.Tg_raster.interp(x=raster.x, y=raster.y)
        _Tg_abs = _Tg_raster.values

        # parameters
        _sza = raster.sza
        _raa = raster.raa
        _vza = raster.vza
        _azi = (180. - _raa) % 360
        _air_mass_ = Misc.air_mass(_sza, _vza)
        _pressure_ = _pressure / self.pressure_ref

        sigma2 = (self.wind_img.wind + 0.586) / 195.3
        _sigma2 = sigma2.interp(x=raster.x, y=raster.y)
        _p_slope_ = self.p_slope(_sza, _vza, _raa, sigma2=_sigma2,
                                 monoview=monoview).values

        # arrays reused at every aot_ref: loaded once
        band = _band_rad.values
        sza_values = _sza.values
        vza_values = _vza.values
        azi_values = _azi.values
        rot_term = rot * np.mean(_pressure_)

        aot_ref_bound = np.min([aot_ref_max_in * 1.02, self.aot_ref_lim])
        aot_refs = np.linspace(0, aot_ref_bound, 21)

        counter = _NegativeCounter(base_threshold=-0.0003)
        xres = []
        for _aot_ref in aot_refs:
            _aot = self.aot_lut.interp(aot_ref=_aot_ref, method='quadratic')  # .values

            arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

            # get LUT values
            _Rdiff = _R_.interp_Rlut_rayleigh(self.szas, sza_values,
                                              self.vzas, vza_values,
                                              self.azis, azi_values,
                                              Nwl, Ny, Nx,
                                              self.Rdiff_lut.sel(aot_ref=_aot_ref, method='nearest').values)

            _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

            #  correction for diffuse light
            Rcorr = band - _Rdiff

            # direct transmittance up/down
            Tdir = np.exp(-(_aot + rot_term) * _air_mass_)  # acutils.Misc.transmittance_dir(_aot, _air_mass_, _rot_raster)

            Rf = np.full((len(iwl_swir), Ny, Nx), np.nan, dtype=np.float32)
            for iwl in iwl_swir:
                if monoview:
                    Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_abs[iwl] * sunglint_eps[iwl] * _p_slope_)
                else:
                    Rf[iwl] = (sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl])
                    Rf[iwl] = Rf[iwl] / (
                            _Tg_abs[iwl] * Tdir[iwl] * sunglint_eps[iwl] * _p_slope_[iwl])

            # Rf[Rf < -0.0004] = np.nan
            # Rf[Rf < 0.] = 0.0
            Rf = np.abs(Rf)
            Rf = np.min(Rf, axis=0)

            Rf = _R_._multiplicate(sunglint_eps, Rf, arr_tmp)
            # slope probability of the reference band: (y, x) for monoview sensors, the same for all the bands
            # (indexing it would take a row of the block)
            _p_slope_ref = _p_slope_ if monoview else _p_slope_[-1]
            Rf = _Tg_abs * Tdir * Rf * _p_slope_ / (sunglint_eps[-1] * _p_slope_ref)

            # sunglint removal
            Rrs_toa = np.asarray((Rcorr - Rf) / np.pi)

            counter.add(Rrs_toa, nwl=7)
            if keep_xres:
                xres.append(xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rrs_toa)),
                                       coords=dict(wl=raster.wl, x=raster.x, y=raster.y, aot_ref=_aot_ref)))

        if keep_xres:
            self.xres = xr.concat(xres, dim='aot_ref')

        aot_ref_maxs = _aot_ref_candidates(aot_refs, counter.counts, neg_pix_max_in, neg_pix_rel_thresh)
        if len(aot_ref_maxs) > 0:
            aot_ref_max = aot_ref_maxs[-1]
        else:
            aot_ref_max = self._cams_aod550.interp(x=raster.x, y=raster.y).mean(['y', 'x']).values

        if update_state:
            self.aot_ref_max = aot_ref_max
        number_of_pixel = raster.water_pixel_number.sum(['y', 'x']).values
        retrieved_aot = xr.Dataset({'aot_ref': (["y", "x"], np.array([[aot_ref_max]])),
                                    'number_of_pixel': (["y", "x"], np.array([[number_of_pixel]]))},
                                   coords=dict(y=[raster.y.mean().values], x=[raster.x.mean().values]))

        return retrieved_aot

    def final_process_chunk(self,
                            raster,
                            aot_ref_raster,
                            iwl_swir=[-2, -1]
                            ):
        r'''
        Rrs, sunglint reflectance and AOT of one block (use with ``xr.map_blocks``, see `final_process`).

        Steps: Rayleigh / aerosol diffuse light correction with the AOT map, sunglint estimated from the
        SWIR bands and removed, conversion from TOA to BOA (division by the total up/down transmittance
        where the reflectance is positive).

        Each pixel uses its own :math:`\tau_{ref}` (interpolated from `aot_ref_raster`) and pressure:

        .. math::

           R_{corr}(\lambda) = R_{TOA}(\lambda) - T_g^{diff}(\lambda)\, \frac{P}{P_0}\, R_{diff}(\lambda;
           \tau_{ref}), \qquad
           T_{dir}(\lambda) = \exp\left[-\left(\tau_a(\lambda) + \tau_R(\lambda)\,\frac{P}{P_0}\right) M\right]

        The sunglint is taken from the first SWIR band :math:`\lambda_{s}` = ``iwl_swir[0]``, with
        :math:`\lambda_{r}` = ``iwl_swir[-1]``, and small negative values are clipped:

        .. math::

           R_f' = \frac{\varepsilon(\lambda_{r})\, P_{slope}(\lambda_{r})\, R_{corr}(\lambda_s)}
           {T_g(\lambda_s)\, T_{dir}(\lambda_s)\, \varepsilon(\lambda_s)\, P_{slope}(\lambda_s)},
           \qquad
           R_f =
           \begin{cases}
              -10^{-5} & \text{if } R_f' < -10^{-3} \\
              0 & \text{if } -10^{-3} \le R_f' < -10^{-5} \\
              R_f' & \text{otherwise}
           \end{cases}

           R_g(\lambda) = T_g(\lambda)\, T_{dir}(\lambda)\, \varepsilon(\lambda)\, P_{slope}(\lambda)\,
           \frac{R_f}{\varepsilon(\lambda_s)\, P_{slope}(\lambda_s)}

        and the remote-sensing reflectance is

        .. math::

           R_{rs}(\lambda) = \frac{R_{corr}(\lambda) - R_g(\lambda)}{\pi\, T_{tot}(\lambda)},
           \qquad
           T_{tot}(\lambda) =
           \begin{cases}
              T_{\downarrow}(\lambda; \theta_s)\, T_{\downarrow}(\lambda; \bar\theta_v)^{1.05}\, T_g(\lambda)
              & \text{if } R_{corr}(\lambda) \ge R_g(\lambda) \\
              1 & \text{otherwise}
           \end{cases}

        ``BRDFg`` is :math:`R_f` and ``aot550`` is :math:`\tau_{ref}`.

        Note: :math:`R_f'` is normalised with :math:`\lambda_r` (last band) but :math:`R_g` with
        :math:`\lambda_s`, so that (without clipping)

        .. math::

           R_g(\lambda_s) = R_{corr}(\lambda_s)\,
           \frac{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)}{\varepsilon(\lambda_s)\, P_{slope}(\lambda_s)}

        which differs from :math:`R_{corr}(\lambda_s)` unless the two products are equal
        (`aerosol_visible_chunk` uses :math:`\lambda_r` for both).

        The wind map ``self.wind_img`` is interpolated to the block; the coarse maps (surface pressure,
        gaseous transmittance, `aot_ref_raster`) are interpolated to the block grid. Nothing is stored on
        the instance.

        :param raster: Dataset block (all wavelengths) with ``bands``, ``sza``, ``raa``, ``vza`` and the
            coordinates ``wl``, ``x``, ``y``; dimensions (wl, y, x)
        :param aot_ref_raster: DataArray (y, x) with the AOT at the reference wavelength, coarse grid,
            in memory
        :param iwl_swir: indices of the two SWIR bands used for the sunglint
        :return: Dataset with ``Rrs`` (wl, y, x), ``BRDFg`` (y, x) and ``aot550`` (y, x), dtype
            ``prod._type``
        '''
        monoview = self.monoview
        _R_ = self._R_
        dtype = self.prod._type

        raster = self._load_block(raster, ('bands', 'sza', 'raa', 'vza'))
        _band_rad = raster.bands
        Nwl, Ny, Nx = _band_rad.shape
        Rrs_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=dtype)
        Rf_tmp = np.full((Ny, Nx), np.nan, dtype=dtype)
        aot_tmp = np.full((Ny, Nx), np.nan, dtype=dtype)

        if Ny > 0 and Nx > 0:
            sunglint_eps = self.sunglint_eps.values
            x_subset, y_subset = _band_rad.x, _band_rad.y
            Ttot_Ed_lut = self.trans_aero_lut.Ttot_Ed

            sigma2 = (self.wind_img.wind + 0.586) / 195.3
            _sigma2_ = sigma2.interp(x=x_subset, y=y_subset)

            arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=np.float32)

            # subsetting
            _sza = raster.sza
            _pressure = self._surface_pressure.interp(x=x_subset, y=y_subset).values
            _Tg_diff_raster = self.Tg_diff_raster.interp(x=x_subset, y=y_subset).values
            _Tg_raster = self.Tg_raster.interp(x=x_subset, y=y_subset).values
            # the AOT map is given at the centres of its coarse pixels: extrapolate it to the pixels of the
            # image edges, outside the first / last centres (NaN otherwise), within the range of the LUT
            _aot_ref_raster = np.clip(aot_ref_raster.interp(x=x_subset, y=y_subset,
                                                            kwargs={'fill_value': 'extrapolate'}).values,
                                      self.aot_refs.min(), self.aot_refs.max())

            _raa = raster.raa
            _vza = raster.vza
            if monoview:
                _vza_mean = _vza.values
            else:
                _vza_mean = np.mean(_vza, axis=0).values

            _azi = (180. - _raa) % 360
            # numpy array: (wl, y, x) for Sentinel-2, (y, x) for monoview sensors (broadcast over wl);
            # a (y, x) DataArray times the (wl, y, x) numpy arrays of Tdir fails with recent xarray
            _air_mass_ = Misc.air_mass(_sza, _vza).values
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

            Rf = np.full((len(iwl_swir), Ny, Nx), np.nan, dtype=np.float32)
            for iwl in iwl_swir:
                if monoview:
                    Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_raster[iwl] * sunglint_eps[iwl] * _p_slope_)
                else:
                    Rf[iwl] = (sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl])
                    Rf[iwl] = Rf[iwl] / (
                            _Tg_raster[iwl] * Tdir[iwl] * sunglint_eps[iwl] * _p_slope_[iwl])

            Rf[Rf < -0.001] = -0.00001  # np.nan
            Rf[Rf < -0.00001] = 0.0
            Rf = Rf[0]  # np.min(Rf, axis=0)
            Rf_tmp[:] = Rf

            Rf = _R_._multiplicate(sunglint_eps, Rf, arr_tmp)
            # slope probability of the reference band: (y, x) for monoview sensors, the same for all the bands
            # (indexing it would take a row of the block)
            _p_slope_ref = _p_slope_ if monoview else _p_slope_[-2]
            Rf = _Tg_raster * Tdir * Rf * _p_slope_ / (sunglint_eps[-2] * _p_slope_ref)

            # sunglint removal
            Rrs_tmp_ = ((Rcorr - Rf) / np.pi)

            # Convert from TOA to BOA for positive values
            Ttot_du[Rrs_tmp_ < 0] = 1.
            Rrs_tmp_ = Rrs_tmp_ / Ttot_du
            Rrs_tmp[:] = Rrs_tmp_
            aot_tmp[:] = _aot_ref_raster

        return xr.Dataset(dict(Rrs=(['wl', "y", "x"], Rrs_tmp),
                               BRDFg=(["y", "x"], Rf_tmp),
                               aot550=(["y", "x"], aot_tmp)
                               ),
                          # plain values: a coordinate DataArray would bring the scalar coordinates of the
                          # block (time, spatial_ref) that the map_blocks template does not have
                          coords=dict(wl=raster.wl.values,
                                      x=raster.x.values,
                                      y=raster.y.values,
                                      ),
                          )

    def final_process(self,
                      raster,
                      aot_ref_raster,
                      iwl_swir=[-2, -1],
                      chunk=512,
                      lazy=True
                      ):
        '''
        Rrs, sunglint reflectance and AOT of the whole image, computed lazily block by block.

        The image is cut into ``chunk x chunk`` blocks (all wavelengths in one block) and
        `final_process_chunk` is mapped on them with ``xr.map_blocks``: the result is a dask-backed
        Dataset, evaluated when it is written or computed, so that neither ``Rrs`` nor any intermediate
        full-size array is in memory. Sets ``self.xres`` (as before) and returns it.

        :param raster: Dataset of the image with ``bands``, ``sza``, ``raa``, ``vza`` (dimensions wl, y, x)
        :param aot_ref_raster: DataArray (y, x) with the AOT at the reference wavelength, on a coarse grid
            (it is loaded in memory if lazy and small)
        :param iwl_swir: indices of the two SWIR bands used for the sunglint
        :param chunk: block size in pixels. A block needs about ``15 * n_bands * chunk**2 * 8`` bytes while it
            is computed (~350 MB for 512, ~1.4 GB for 1024, per thread). 512 divides the 1024 shards of
            ``export_to_zarr``, which merges four blocks per shard without recomputing anything.
        :param lazy: if False, compute the result immediately (previous behaviour, needs the memory of the
            full output)
        :return: Dataset ``xres`` with ``Rrs`` (wl, y, x), ``BRDFg`` (y, x), ``aot550`` (y, x)

        Notes
        -----
        The blocks are computed later, when the result is written, so the function works on a frozen copy of
        the kernel taken now (`_snapshot`): changing ``kernel.wind_img`` or recomputing the LUTs after this
        call does not change the result.
        '''
        variables = [v for v in ('bands', 'sza', 'raa', 'vza') if v in raster.data_vars]
        block = raster[variables].transpose('wl', 'y', 'x').chunk({'y': chunk, 'x': chunk, 'wl': -1})
        aot_ref_raster = _load_small(aot_ref_raster)

        # structure of the result: same blocks as the bands, no extra coordinate
        dtype = self.prod._type
        bands = block['bands'].reset_coords(drop=True)
        plane = bands.isel(wl=0, drop=True)
        template = xr.Dataset({'Rrs': xr.full_like(bands, np.nan, dtype=dtype),
                               'BRDFg': xr.full_like(plane, np.nan, dtype=dtype),
                               'aot550': xr.full_like(plane, np.nan, dtype=dtype)})

        xres = xr.map_blocks(self._snapshot().final_process_chunk, block,
                             kwargs=dict(aot_ref_raster=aot_ref_raster, iwl_swir=iwl_swir),
                             template=template)
        self.xres = xres if lazy else xres.compute()
        return self.xres

    def _snapshot(self):
        '''
        Frozen copy of what the block functions read, without the heavy inputs.

        A lazy result is computed after the call that built it, when the kernel may have changed (e.g. the
        wind map is filtered, the LUTs are prepared again). The copy keeps the attributes as they are now
        (they are re-assigned, not modified in place, by the kernel; ``wind_img`` is copied) and drops the
        raster, the raw LUTs and the CAMS data, so the graph does not carry them (nor pickle them).

        :return: shallow copy of the kernel, usable for ``final_process_chunk`` only
        '''
        state = copy.copy(self)
        _ = self._surface_pressure  # load now: the copy has no access to the CAMS data
        state._cache = dict(self._cache)
        state.wind_img = self.wind_img.copy(deep=False)
        state.prod = types.SimpleNamespace(_type=self.prod._type, sensor=self.prod.sensor)
        for name in ('raster', 'aero_lut', 'trans_lut', 'cams', 'coarse_masked_raster',
                     'xres', 'aot_ref_img', 'aot_ref_raster'):
            state.__dict__.pop(name, None)
        return state

    @staticmethod
    @njit
    def filter2d(image, weight, windows):
        r'''
        Weighted mean filter of an image (numba).

        Each pixel with a positive weight is replaced by the weighted mean of its neighbours in a
        ``windows`` window (pixels with a null / NaN weight are ignored):

        .. math::

           \tilde X(i, j) = \frac{\sum_{(k, l) \in \mathcal{W}_{ij},\, w_{kl} > 0} w_{kl}\, X(k, l)}
           {\sum_{(k, l) \in \mathcal{W}_{ij},\, w_{kl} > 0} w_{kl}}

        The window is truncated at the edges of the image; pixels with a null / NaN weight keep their value.
        All the means are computed from the input values (`image` is not modified).

        :param image: 2D array to smooth
        :param weight: 2D array of weights (e.g. inverse uncertainty), same shape
        :param windows: size of the window (rows, columns)
        :return: smoothed copy of `image`
        '''
        M, N = np.shape(image)
        Mf, Nf = windows
        Mf2 = Mf // 2
        Nf2 = Nf // 2
        threshold = 0
        result = image.copy()
        for i in range(M):
            for j in range(N):
                num = 0.0
                norm = 0.0
                if weight[i, j] > threshold:
                    for ii in range(Mf):
                        ix = i - Mf2 + ii
                        if 0 <= ix < M:
                            for jj in range(Nf):

                                iy = j - Nf2 + jj
                                if 0 <= iy < N:
                                    wgt = weight[ix, iy]
                                    if wgt > 0.:
                                        num += (wgt * image[ix, iy])
                                        norm += wgt
                    result[i, j] = num / norm
        return result

    @staticmethod
    def conv_mapping(x):
        """
        Nan-mean convolution: function for ``scipy.ndimage.generic_filter``.

        The central value is replaced by the mean of the other values if it is NaN, kept if all the others
        are NaN, and else averaged with them.

        :param x: flattened values of the footprint (the central one is ``x[len(x) // 2]``)
        :return: filtered value
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
        '''
        Smooth one variable of a (small) raster, e.g. the wind map: weighted mean (`filter2d`, weights
        ``1`` where the variable is finite and non null) then NaN-aware mean (`conv_mapping`).

        The raster is loaded in memory: use it on the coarse wind map, not on a full-resolution image.

        :param raster: Dataset (y, x) containing `varname`
        :param varname: name of the variable
        :param windows: window of `filter2d` (rows, columns)
        :param mask: footprint of the second filter
        :return: Dataset with the smoothed variable ``varname`` and the coordinates ``y``, ``x``
        '''
        weights = 1 + 0 * (1 / raster[varname] ** 2).values
        param = raster[varname].values
        raster_smoothed = self.filter2d(param, weights, windows)
        if np.all(mask == 1):
            res = self._nanmean_box(raster_smoothed, mask.shape)
        else:
            res = ndimage.generic_filter(raster_smoothed, function=self.conv_mapping, footprint=mask,
                                         mode='nearest')

        return xr.Dataset({varname: (["y", "x"], res)}, coords=dict(y=raster.y, x=raster.x))

    @staticmethod
    def _nanmean_box(image, size):
        '''
        NaN-aware mean over a rectangular window, edges extended by the nearest value.

        Same result as ``generic_filter(image, conv_mapping, footprint=np.ones(size), mode='nearest')``
        (to floating point rounding: the central value is NaN or not, the mean of the valid values of the
        window is returned, NaN if there is none), but computed with two box filters instead of a Python call
        per pixel and per window element: for a 256 x 256 window it is thousands of times faster.

        :param image: 2D array
        :param size: window shape (rows, columns)
        :return: filtered array, float64
        '''
        valid = ~np.isnan(image)
        total = ndimage.uniform_filter(np.where(valid, image, 0.), size=size, mode='nearest')
        count = ndimage.uniform_filter(valid.astype(np.float64), size=size, mode='nearest')
        out = np.full(image.shape, np.nan)
        ok = count > 0.5 / (size[0] * size[1])  # a true count is at least 1 / (window size), the rest is rounding
        out[ok] = total[ok] / count[ok]
        return out

    # ------------------------------------------------------------------
    # LEGACY eager versions, kept unchanged for backward compatibility (not dask-friendly).
    # They allocate full-size numpy arrays and one full raster per aot_ref, and write results to self.
    # Superseded by: rugosity_est_chunk, aerosol_swir_chunk, aerosol_visible_chunk, final_process.
    # `rugosity_est` is known not to work (see its TODO). Delete them once nothing calls them.
    # ------------------------------------------------------------------

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
            self.aot_ref_max = aot_refs[-1]  # self.aot_ref_cams_max

        # for further lut interpolation:
        # self.aot_ref_max = np.max([0.011, aot_ref_max])

        # to stay within LUT aot range
        if ~((self.aot_ref_max > 0) & (self.aot_ref_max < self.aot_ref_lim)):
            self.aot_ref_max = np.min([self.aot_ref_cams_max, self.aot_ref_lim])

        # TODO uncomment after end of dev
        # del self.xres

    def aerosol_visible_chunk_old(self,
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

        aot_ref_bound = np.min([self.aot_ref_max * 1.02, self.aot_ref_lim])
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

            Rf = np.full((len(iwl_swir), Ny, Nx), np.nan, dtype=np.float32)
            for iwl in iwl_swir:
                if monoview:
                    Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_abs[iwl] * sunglint_eps[iwl] * _p_slope_)
                else:
                    Rf[iwl] = (sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl])
                    Rf[iwl] = Rf[iwl] / (
                            _Tg_abs[iwl] * Tdir[iwl] * sunglint_eps[iwl] * _p_slope_[iwl])

            #Rf[Rf < -0.0004] = np.nan
            #Rf[Rf < 0.] = 0.0
            Rf =np.abs(Rf)
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
        metric = metric.where(metric.isel(aot_ref=0) > -0.0003)

        metric_pixnum = metric.where(metric < 0).count(['x', 'y'])
        #print(metric_pixnum)
        neg_pix_max = np.max([self.neg_pix_max, metric_pixnum.max()])
        metric_pixnum_max = neg_pix_max * neg_pix_rel_thresh

        aot_ref_maxs = metric_pixnum.where(metric_pixnum < metric_pixnum_max, drop=True).aot_ref
        if len(aot_ref_maxs) > 0:
            aot_ref_max = aot_ref_maxs.values[-1]
        else:
            retrieved_aot = (self.cams.cams_aod.interp(wl=550, method='quadratic')).interp(x=raster.x, y=raster.y)
            retrieved_aot.name = 'aot_ref'
            return retrieved_aot

        self.aot_ref_max = aot_ref_max
        print('----------------')
        print(aot_ref_max)
        print('----------------')

        # ind = np.nanargmin((aot_refs - aot_ref_max) ** 2)
        # map_pixels = metric.sel(aot_ref=aot_refs[ind])

        map_pixels = self.cams.cams_aod.interp(wl=550, method='quadratic'
                                               ).interp(x=metric.x, y=metric.y)

        map_pixels = map_pixels.where(map_pixels < 0)

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