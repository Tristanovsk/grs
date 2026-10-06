'''
Pixel classification (cloud, water, cirrus, ... masks) for large rasters, computed in parallel with dask.

Overview
========

:class:`Masking` builds a single ``uint32`` bit-field raster, ``flags``, in which every bit is
one binary test on the image (see the table below), then derives a binary good/bad ``mask``
from selected bits with :meth:`Masking.create_mask`.

Typical use::

    prod = Product('scene_L2A.nc')                      # dask-backed, chunks (wl=-1, y=1024, x=1024)
    masking = Masking(prod.raster, tile_size=1024)      # keep tile_size equal to the y/x chunk size
    raster = masking.process(output='prod')             # adds the variable `flags` and `flag_*` statistics
    mask = Masking.create_mask(raster.flags,
                               tomask=[0, 1, 2, 7],     # bits that make a pixel bad
                               tokeep=[3])              # pixel must have (at least one of) these bits set
    raster['mask'] = mask

Flag bits
=========

====  ===========================  =============================================================
bit   name                         meaning
====  ===========================  =============================================================
0     ``nodata``                   NaN or crazy value (> 1e3) in the input
1     ``cloud_p06``                s2cloudless, low confidence (probability > 0.6, dilation 12)  [Sentinel-2]
2     ``cloud_p08``                s2cloudless, high confidence (probability > 0.8, dilation 5)  [Sentinel-2]
3     ``water_swir_visible_index`` (visible - SWIR) / (visible + SWIR) above threshold
4     ``water_red_visible_index``  (visible - NIR) / (visible + NIR) above threshold (NDWI-like)
5     ``thin_cirrus``              cirrus band above the thin cirrus threshold
6     ``opac_cirrus``              cirrus band above the opaque cirrus threshold
7     ``high_swir``                SWIR2 reflectance above threshold (bright cloud, glint, ...)
8     ``surfwater_land``           input variable ``surfwater`` == 0
9     ``surfwater_water``          input variable ``surfwater`` == 1
10    ``surfwater_cloud_and_shadow``  input variable ``surfwater`` >= 2
11    ``OCM_Thick_Cloud``
12
13
14-20  ``l1_*``                     Landsat level-1 QA bits copied from ``flags_l1``  [Landsat only]
====  ===========================  =============================================================

How it is computed (dask)
=========================

* The flags are built as **one lazy dask graph**: the band cube is cut into ``tile_size`` x
  ``tile_size`` tiles (all wavelengths in one chunk) and every tile goes through all the masks in
  a single task. Tiles are processed in parallel by the active dask scheduler (threads by default,
  or a ``dask.distributed`` cluster if a client exists); there is no Python tile loop.
* Only the bands a mask needs are read (band selection is fused into the read when the dataset
  chunks match ``tile_size``).
* Neighborhood operations (s2cloudless smoothing and dilation) use overlapping tiles
  (``dask.array.overlap``) with a halo larger than their radius, so the result is identical to a
  run on the full image while the probability maps only ever exist for one tile.
* :meth:`Masking.process` evaluates flags **and** per-flag statistics in one ``dask.compute``
  call, so every tile is processed exactly once. The result (``uint32``, 4 bytes per pixel) is kept
  in memory and stored in the dataset: cloud detection is never repeated downstream.
* Peak memory is about ``n_workers x tile_size**2 x (n_bands x 4 bytes x ~3)``. Lower
  ``tile_size`` or the number of workers (``dask.config.set(num_workers=...)``) if memory is tight.
  s2cloudless uses LightGBM, which has its own threads: on many cores set ``OMP_NUM_THREADS`` or use
  a distributed cluster with one thread per worker to avoid oversubscription.
'''

import logging
from collections import namedtuple
from functools import lru_cache

import dask
import dask.array as da
from dask.array.overlap import overlap, trim_internal
import numpy as np
import pandas as pd
import xarray as xr

import omnicloudmask as OCM
from s2cloudless import S2PixelCloudDetector

FLAG_NAME = 'flags'


class Settings:
    '''
    Settings used for the different masking procedures.

    Attributes
    ----------
    bvis, bnir, bswir, bswir2 : float
        Central wavelengths (nm) of the visible (490), NIR (842), SWIR (1610) and SWIR2 (2190) bands.
        The nearest available band is used, except for ``bswir2`` which must exist exactly.
    bcirrus : float or None
        Cirrus band (nm); set from the constellation (1375 for Sentinel-2, 1370 for Landsat 8/9).
    bwv : float
        Water-vapour band (nm); currently not used by any mask.
    swir_threshold : float
        ``high_swir`` flag: SWIR2 reflectance above this value.
    ndwi_threshold : float
        ``water_red_visible_index`` flag: (visible - NIR) / (visible + NIR) above this value.
    vis_swir_index_threshold : float
        ``water_swir_visible_index`` flag: (visible - SWIR) / (visible + SWIR) above this value.
    thin_cirrus_threshold, opac_cirrus_threshold : float
        Cirrus-band reflectance thresholds of the thin and opaque cirrus flags.
    cloud_average_over : int
        s2cloudless: size of the averaging window applied to the probability maps (pixels).
    low_confid_cloud_proba_thresh, high_confid_cloud_proba_thresh : float
        s2cloudless probability thresholds of the low and high confidence cloud flags.
    low_confid_cloud_dilation, high_confid_cloud_dilation : int
        s2cloudless dilation radius (pixels) of the low and high confidence cloud flags.
    number_of_flags : int
        Size of the flag name/description tables (bits of the ``uint32`` flags).
    '''

    def __init__(self):
        self.bvis = 490
        self.bnir = 842
        self.bswir = 1610
        self.bswir2 = 2190

        # set extra bands (e.g., cirrus, water vapor)
        self.bcirrus = None
        self.bwv = 945

        # mask thresholding parameters
        self.swir_threshold = 0.2
        self.ndwi_threshold = 0.0
        self.vis_swir_index_threshold = 0.2
        self.thin_cirrus_threshold = 3e-3
        self.opac_cirrus_threshold = 9e-3

        # settings for s2cloudless masks
        self.cloud_average_over = 1
        self.low_confid_cloud_proba_thresh = 0.6
        self.low_confid_cloud_dilation = 12
        self.high_confid_cloud_proba_thresh = 0.8
        self.high_confid_cloud_dilation = 5

        # settings OmniCloudMask
        self.OCM_Rband = 665
        self.OCM_Gband = 560
        self.OCM_NIRband = 865
        self.OCM_batch_size = 1
        self.OCM_inference_dtype = "fp32"
        self.OCM_inference_device = 'cpu'
        self.OCM_mosaic_device = 'cpu'

        # general
        self.number_of_flags = 32


# One masking test, applied to every tile:
#   fn(tile)  ORs its bits into tile.f
#   bands     indices (along 'wl') of the bands it reads; None = all bands, [] = none
#   halo      extra pixels around the tile it needs (neighbourhood operations)
#   aux       names of the 2-D dataset variables it reads (e.g. 'surfwater')
_Step = namedtuple('_Step', 'fn bands halo aux')


class _Tile:
    '''Numpy data handed to a step for one (possibly halo-padded) tile.

    ``cube`` is the band cube (wl, y, x) restricted to the needed bands (or None), ``aux`` maps
    variable names to 2-D arrays, ``f`` is the ``uint32`` flag array of the tile, all with the
    same (y, x) shape.
    '''
    __slots__ = ('cube', 'slot', 'aux', 'f')

    def __init__(self, cube, slot, aux, f):
        self.cube, self.slot, self.aux, self.f = cube, slot, aux, f

    def band(self, iwl):
        '''View of the band with index `iwl` (index in the full wl axis of the dataset).'''
        return self.cube[iwl if self.slot is None else self.slot[iwl]]


def _set(flags, bit, cond):
    '''In-place: set bit `bit` of `flags` where `cond` is True (no temporary array).'''
    np.bitwise_or(flags, np.uint32(1 << bit), out=flags, where=cond)


def _chunks(n, size, minimum=0):
    '''Regular chunk sizes along an axis of length `n`; a last chunk smaller than `minimum` is merged
    into the previous one (dask overlaps need chunks at least as large as the halo).'''
    sizes = [size] * (n // size) + ([n % size] if n % size else [])
    if len(sizes) > 1 and sizes[-1] < max(minimum, 1):
        sizes[-2] += sizes.pop()
    return tuple(sizes)


@lru_cache(maxsize=8)
def _s2_detector(threshold, average_over, dilation_size):
    '''s2cloudless detector, created once per process (workers) instead of being shipped with every task.'''
    return S2PixelCloudDetector(threshold=threshold, average_over=average_over,
                                dilation_size=dilation_size, all_bands=True)


def _make_block_fn(steps, slot, has_cube, aux_names):
    '''Build the function applied by dask to each tile: runs all `steps` and returns the tile flags.'''

    def block_fn(*arrays):
        arrays = list(arrays)
        cube = arrays.pop(0) if has_cube else None
        aux = dict(zip(aux_names, arrays))
        shape = cube.shape[1:] if has_cube else arrays[0].shape
        tile = _Tile(cube, slot, aux, np.zeros(shape, dtype=np.uint32))
        for step in steps:
            step.fn(tile)
        return tile.f

    return block_fn


class Masking(Settings):
    '''
    Pixel classification and masking of an image raster.

    Computes the bit-field ``flags`` (see the module documentation for the bit table) and the
    derived binary mask. The flags are built with dask, tile by tile, in parallel.

    Parameters
    ----------
    prod : xarray.Dataset
        Image dataset. Required: variable ``bands`` with dimensions (wl, y, x) in any order and
        the attribute ``constellation_id`` (``'S2...'``, ``'L8'`` or ``'L9'``). Optional variables:
        ``surfwater`` (2-D) for the surfwater flags, ``flags_l1`` (2-D integer) for Landsat.
        The dataset may be numpy- or dask-backed (for dask, use chunks ``wl=-1`` and ``y``/``x`` equal
        to `tile_size`, so no rechunking is needed).
    tile_size : int, optional
        Edge (pixels) of the tiles processed by one task (default 1024). It must be at least as large
        as the s2cloudless halo (``max(dilation) + average_over + 1`` pixels).

    Attributes
    ----------
    flags : numpy.ndarray or None
        ``uint32`` flags (y, x), available after `process` or any single-mask method.
    flag_names, flag_descriptions : numpy.ndarray of object
        Name and description of every bit (None for unused bits).
    flag_stats : dict
        Fraction of flagged pixels per flag, keys ``'flag_<name>'``.
    sentinel2 : bool
        True for Sentinel-2 (s2cloudless is used), False for Landsat (level-1 QA bits are copied).
    '''

    def __init__(self, prod: xr.Dataset, tile_size: int = 1024):
        Settings.__init__(self)
        self.prod = prod
        y_size = prod.sizes["y"]
        x_size = prod.sizes["x"]
        self.tile_size = int(tile_size)

        self.sentinel2 = False
        if "S2" in prod.attrs['constellation_id']:
            self.bcirrus = 1375
            self.sentinel2 = True
        elif (prod.attrs['constellation_id'] == 'L8') or (prod.attrs['constellation_id'] == 'L9'):
            self.bcirrus = 1370

        self.flag_descriptions = np.empty(self.number_of_flags, dtype='object')
        self.flag_names = np.empty(self.number_of_flags, dtype='object')

        # Allocate by default in memory as uint32
        # np.zeros fills it with 0 (clear sky default), or use np.empty
        self.flags = np.zeros((y_size, x_size), dtype=np.uint32)
        self.flag_stats = {}

    # ------------------------------------------------------------------
    # machinery
    # ------------------------------------------------------------------
    def _register(self, bit, name, description):
        '''Record the name and description of flag bit `bit`.'''
        self.flag_names[bit] = name
        self.flag_descriptions[bit] = description

    def _wl_index(self, wl, nearest=True):
        '''Position of wavelength `wl` along the wl axis (nearest, or exact like ``sel(wl=wl)``).'''
        index = self.prod.indexes['wl']
        if nearest:
            return int(index.get_indexer([wl], method='nearest')[0])
        return int(index.get_loc(wl))

    def _require(self, name):
        '''Return dataset variable `name`, failing early with a clear message if it is missing.'''
        if name not in self.prod:
            raise KeyError(f"the dataset has no variable '{name}', needed by this mask")
        return self.prod[name]

    def _flags_lazy(self, steps):
        '''Lazy dask array (y, x, uint32) with the flags of `steps`.

        All steps run in one task per tile. Tiles overlap by the largest halo of the steps, are
        processed whole (halo included) and trimmed afterwards.
        '''
        ny, nx = self.prod.sizes['y'], self.prod.sizes['x']
        halo = max(s.halo for s in steps)
        if self.tile_size < halo:
            raise ValueError(f'tile_size ({self.tile_size}) must be at least the halo ({halo} pixels)')
        chunks = {'y': _chunks(ny, self.tile_size, halo), 'x': _chunks(nx, self.tile_size, halo)}

        # bands needed by the steps
        if any(s.bands is None for s in steps):
            needed, slot = None, None  # all bands
        else:
            needed = sorted({i for s in steps for i in s.bands})
            slot = {iwl: k for k, iwl in enumerate(needed)}
        has_cube = needed != []
        aux_names = sorted({n for s in steps for n in s.aux})

        arrays, depths = [], []
        if has_cube:
            bands = self._require('bands').transpose('wl', 'y', 'x').chunk({'wl': -1, **chunks})
            if needed is not None:
                bands = bands.isel(wl=needed)  # fused into the read: other bands are not loaded
            arrays.append(bands.data)
            depths.append({0: 0, 1: halo, 2: halo})
        for name in aux_names:
            arrays.append(self._require(name).transpose('y', 'x').chunk(chunks).data)
            depths.append({0: halo, 1: halo})

        if halo > 0:
            arrays = [overlap(a, depth=d, boundary='none') for a, d in zip(arrays, depths)]

        flags = da.map_blocks(_make_block_fn(steps, slot, has_cube, aux_names), *arrays,
                              drop_axis=0 if has_cube else [],
                              dtype=np.uint32,
                              meta=np.empty((0, 0), dtype=np.uint32))
        if halo > 0:
            flags = trim_internal(flags, {0: halo, 1: halo}, boundary='none')
        return flags

    def _run(self, steps):
        '''Compute `steps` (in parallel) and OR their bits into ``self.flags``.'''
        result = self._flags_lazy(steps).compute()
        if self.flags is None:
            self.flags = result
        else:
            np.bitwise_or(self.flags, result, out=self.flags)

    # ------------------------------------------------------------------
    # masks: *_step builds the step (and registers the flag names),
    #        the public method runs it alone on the whole image
    # ------------------------------------------------------------------
    def _nodata_step(self, bitmask=0, name='nodata', description='nodata in input image'):
        '''Step: flag nodata, i.e. NaN or crazy values (> 1e3) in the band with index 1.'''
        self._register(bitmask, name, description)

        def fn(t):
            b = t.band(1)
            _set(t.f, bitmask, np.isnan(b) | (b > 1e3))

        return _Step(fn, [1], 0, ())

    def _s2cloudless_step(self,
                          bitmasks=(1, 2),
                          names=('cloud_p06', 'cloud_p08'),
                          descriptions=('low confidence cloud mask from s2cloudless with settings proba.',
                                        'high confidence cloud mask from s2cloudless with settings proba.')):
        '''Step: s2cloudless cloud masks with two levels of confidence (needs all bands).

        See https://github.com/sentinel-hub/sentinel2-cloud-detector. The probability map is computed
        once and thresholded/dilated twice. The tiles carry a halo larger than the smoothing and
        dilation radii, so the result matches a full-image run.
        '''
        low_args = (self.low_confid_cloud_proba_thresh, self.cloud_average_over,
                    self.low_confid_cloud_dilation)
        high_args = (self.high_confid_cloud_proba_thresh, self.cloud_average_over,
                     self.high_confid_cloud_dilation)
        self._register(bitmasks[0], names[0], descriptions[0] + ' threshold={:.2f}, dilation_size={:d}'.format(
            self.low_confid_cloud_proba_thresh, self.low_confid_cloud_dilation))
        self._register(bitmasks[1], names[1], descriptions[1] + ' threshold={:.2f}, dilation_size={:d}'.format(
            self.high_confid_cloud_proba_thresh, self.high_confid_cloud_dilation))
        halo = max(self.low_confid_cloud_dilation, self.high_confid_cloud_dilation) + self.cloud_average_over + 1

        def fn(t):
            low, high = _s2_detector(*low_args), _s2_detector(*high_args)
            data = np.ascontiguousarray(np.moveaxis(t.cube, 0, -1))[np.newaxis]  # (1, h, w, wl)
            prob = low.get_cloud_probability_maps(data)
            del data
            for detector, bit in ((low, bitmasks[0]), (high, bitmasks[1])):
                _set(t.f, bit, detector.get_mask_from_prob(prob)[0] == 1)

        return _Step(fn, None, halo, ())

    # TODO process downscaled image to save time (e.g. resolution 60 m)
    def _get_omnicloudmask(self,
                      bitmasks=(11, 12, 13),
                      names=('OCM_Thick_Cloud', 'OCM_Thin_Cloud', 'OCM_Cloud_Shadow'),
                      descriptions=('Thick Cloud from OmniCloudMask version ' + OCM.__version__,
                                    'Thin Cloud from OmniCloudMask version ' + OCM.__version__,
                                    'Cloud Shadow from OmniCloudMask version ' + OCM.__version__)):
        '''Apply OmniCloudMAsk for clouds and cloud shadows masking

        Outputs:
            0 = Clear
            1 = Thick Cloud
            2 = Thin Cloud
            3 = Cloud Shadow

        see https://github.com/DPIRD-DMA/OmniCloudMask

        refs:
         Wright, N., Duncan, J. M. A., Callow, J. N., Thompson, S. E., & George, R. J. (2025).
         Training sensor-agnostic deep learning models for remote sensing:
         Achieving state-of-the-art cloud and cloud shadow identification with OmniCloudMask.
         Remote Sensing of Environment, 322, 114694. https://doi.org/10.1016/J.RSE.2025.114694

        :param rgnir: raster xarray object with the red, green and nir bands
        :return omnimask: raster of the retrieved mask
        '''
        rgnir = self.prod.bands.sel(wl=[self.OCM_Rband, self.OCM_Gband, self.OCM_NIRband])
        self.rgnir = rgnir
        for ii in range(len(bitmasks)):
            self._register(bitmasks[ii], names[ii], descriptions[ii])

        logging.info('OmniCloudMask processing')
        omnimask = OCM.predict_from_array(rgnir.fillna(0).values,
                                          batch_size=self.OCM_batch_size,
                                          inference_dtype=self.OCM_inference_dtype,
                                          inference_device=self.OCM_inference_device,
                                          mosaic_device=self.OCM_mosaic_device
                                          )[0]
        self.omnimask=omnimask
        # set flags raster
        self.flags = self.flags + (
                ((1 == omnimask) << bitmasks[0]) +
                ((2 == omnimask) << bitmasks[1]) +
                ((3 == omnimask) << bitmasks[2])
        )

    def _water_step(self,
                    bitmasks=(3, 4),
                    names=('water_swir_visible_index', 'water_red_visible_index'),
                    descriptions=('water mask from normalized index from swir and visible band',
                                  'water mask from normalized index from nir and visible band, '
                                  'warning could fail for turbid waters')):
        '''Step: water/land masks from the normalised indices (visible, SWIR) and (visible, NIR).'''
        ivis, inir, iswir = (self._wl_index(w) for w in (self.bvis, self.bnir, self.bswir))
        self._register(bitmasks[0], names[0], descriptions[0] +
                       ', bands centered on {:.1f} and {:.1f} nm'.format(self.bvis, self.bswir) +
                       ', threshold={:.2f}'.format(self.vis_swir_index_threshold))
        self._register(bitmasks[1], names[1], descriptions[1] +
                       ', bands centered on {:.1f} and {:.1f} nm'.format(self.bvis, self.bnir) +
                       ', threshold={:.2f}'.format(self.ndwi_threshold))

        def fn(t):
            vis = t.band(ivis)
            with np.errstate(divide='ignore', invalid='ignore'):
                other = t.band(iswir)
                _set(t.f, bitmasks[0], (vis - other) / (vis + other) > self.vis_swir_index_threshold)
                other = t.band(inir)
                _set(t.f, bitmasks[1], (vis - other) / (vis + other) > self.ndwi_threshold)

        return _Step(fn, [ivis, inir, iswir], 0, ())

    def _cirrus_step(self,
                     bitmasks=(5, 6),
                     names=('thin_cirrus', 'opac_cirrus'),
                     descriptions=('thin cirrus mask from cirrus band',
                                   'opac cirrus mask from cirrus band')):
        '''Step: thin and opaque high clouds from the spectral band around 1375 nm.'''
        icirrus = self._wl_index(self.bcirrus)
        thresholds = (self.thin_cirrus_threshold, self.opac_cirrus_threshold)
        for ii in range(2):
            self._register(bitmasks[ii], names[ii], descriptions[ii] + ', threshold {:.4f}'.format(thresholds[ii]))

        def fn(t):
            cirrus = t.band(icirrus)
            for ii in range(2):
                _set(t.f, bitmasks[ii], cirrus > thresholds[ii])

        return _Step(fn, [icirrus], 0, ())

    def _high_swir_step(self,
                        bitmask=7,
                        name='high_swir',
                        description='high swir (bright cloud, too bright reflection...)'):
        '''Step: (too) high values of the SWIR2 band, to remove potential cloud contamination.'''
        iswir2 = self._wl_index(self.bswir2, nearest=False)
        self._register(bitmask, name, description + ', threshold {:.4f}'.format(self.swir_threshold))

        def fn(t):
            _set(t.f, bitmask, t.band(iswir2) > self.swir_threshold)

        return _Step(fn, [iswir2], 0, ())

    def _surfwater_step(self,
                        bitmasks=(8, 9, 10),
                        names=('surfwater_land', 'surfwater_water', 'surfwater_cloud_and_shadow'),
                        descriptions=('land mask from surfwater input file',
                                      'water mask from surfwater input file',
                                      'cloud and shadow mask from surfwater input file')):
        '''Step: land / water / cloud-and-shadow classes of the 2-D input variable ``surfwater``.'''
        for ii in range(3):
            self._register(bitmasks[ii], names[ii], descriptions[ii])
        self._require('surfwater')

        def fn(t):
            sw = t.aux['surfwater']
            _set(t.f, bitmasks[0], sw == 0)
            _set(t.f, bitmasks[1], sw == 1)
            _set(t.f, bitmasks[2], sw >= 2)

        return _Step(fn, [], 0, ('surfwater',))

    def _landsat_step(self,
                      bitmasks=(14, 15, 16, 17, 18, 19, 20),
                      landsat_bitmasks=(1, 2, 3, 4, 5, 6, 7),
                      names=('l1_dilated_cloud', 'l1_cirrus', 'l1_cloud',
                             'l1_cloud_shadow', 'l1_snow', 'l1_clear', 'l1_water'),
                      descriptions=('landsat dilated cloud mask',
                                    'landsat cirrus mask',
                                    'landsat cloud mask',
                                    'landsat cloud shadow mask',
                                    'landsat snow mask',
                                    'landsat clear mask',
                                    'landsat water mask'),
                      flags_landsat_band='flags_l1'):
        '''Step: copy the main Landsat level-1 QA bits (integer variable `flags_landsat_band`) into the flags.'''
        for ii, bit in enumerate(bitmasks):
            self._register(bit, names[ii], descriptions[ii])
        self._require(flags_landsat_band)

        def fn(t):
            f1 = t.aux[flags_landsat_band]
            for ii, bit in enumerate(bitmasks):
                _set(t.f, bit, (f1 & (1 << landsat_bitmasks[ii])) != 0)

        return _Step(fn, [], 0, (flags_landsat_band,))

    # public single-mask API (same names as before; each one runs its own pass over the image
    # and ORs its bits into `self.flags`; prefer `process()`, which does everything in one pass)
    def nodata_mask(self, *args, **kwargs):
        '''Compute the nodata flag (bit 0) and add it to ``self.flags``.'''
        self._run([self._nodata_step(*args, **kwargs)])

    def s2cloudless_mask(self, *args, **kwargs):
        '''Compute the s2cloudless cloud flags (bits 1, 2) and add them to ``self.flags``.'''
        logging.info('cloud masking with s2cloudless')
        self._run([self._s2cloudless_step(*args, **kwargs)])

    def water_mask(self, *args, **kwargs):
        '''Compute the water index flags (bits 3, 4) and add them to ``self.flags``.'''
        logging.info('water masking')
        self._run([self._water_step(*args, **kwargs)])

    def cirrus_mask(self, *args, **kwargs):
        '''Compute the cirrus flags (bits 5, 6) and add them to ``self.flags``.'''
        logging.info('cirrus masking')
        self._run([self._cirrus_step(*args, **kwargs)])

    def high_swir_mask(self, *args, **kwargs):
        '''Compute the high SWIR flag (bit 7) and add it to ``self.flags``.'''
        logging.info('high swir masking')
        self._run([self._high_swir_step(*args, **kwargs)])

    def surfwater_mask(self, *args, **kwargs):
        '''Compute the surfwater flags (bits 8-10) and add them to ``self.flags``.'''
        logging.info('surfwater masking')
        self._run([self._surfwater_step(*args, **kwargs)])

    def duplicate_landsat_flags(self, *args, **kwargs):
        '''Copy the Landsat level-1 QA bits (bits 11-17) into ``self.flags``.'''
        logging.info('duplicating landsat level 1 flags')
        self._run([self._landsat_step(*args, **kwargs)])

    # ------------------------------------------------------------------
    # statistics
    # ------------------------------------------------------------------
    def _used_bits(self):
        '''Bits that have a registered name.'''
        return [ii for ii, name in enumerate(self.flag_names) if name is not None and name != 'None']

    def _store_stats(self, bits, counts, size):
        '''Save the fraction of flagged pixels of each bit in `flag_stats` and in the dataset attributes.'''
        for bit, count in zip(bits, counts):
            self.flag_stats['flag_' + str(self.flag_names[bit])] = float(count) / size
        self.prod.attrs.update(self.flag_stats)

    def get_stats(self):
        '''
        Compute image statistics for each flag (fraction of flagged pixels) from ``self.flags`` and save
        them in ``flag_stats`` and in ``prod.attrs``. Counted block-wise: no full-size temporary.
        (`process` obtains the same numbers during its single dask pass; call this after the
        single-mask methods.)
        '''
        flags = self.flags
        bits = self._used_bits()
        counts = np.zeros(len(bits), dtype=np.int64)
        rows = max(1, (1 << 22) // max(flags.shape[1], 1))
        for y0 in range(0, flags.shape[0], rows):
            block = flags[y0:y0 + rows]
            for k, bit in enumerate(bits):
                counts[k] += np.count_nonzero(block & np.uint32(1 << bit))
        self._store_stats(bits, counts, flags.size)

    def print_stats(self):
        '''
        Provide pandas Dataframe with image statistics of each flag.

        Returns
        -------
        pandas.DataFrame
            Columns ``name``, ``description``, ``bit`` and ``statistics`` (fraction of flagged pixels).
        '''

        pflags = self.prod[FLAG_NAME]

        # construct dataframe:
        dflags = pd.DataFrame({'name': pflags.attrs["flag_names"]})
        dflags['description'] = pflags.attrs["flag_descriptions"]
        dflags['bit'] = dflags.index
        dflags = dflags[dflags.name != "None"]
        stats = []
        for name in dflags.name:
            stats.append(self.prod.attrs["flag_" + name])
        dflags['statistics'] = stats
        return dflags

    def process(self,
                output="prod"):
        '''
        Generate the flags raster and attributes in ONE pass over the tiles.

        Builds the lazy dask graph of all masks (nodata, cloud or Landsat QA, water, cirrus if a cirrus
        band exists, high SWIR, surfwater), then evaluates the flags and the per-flag statistics together
        with a single ``dask.compute``: every tile is processed once, in parallel. Previously computed
        flags are discarded. The flags are stored (numpy-backed) in ``self.flags`` and in the dataset
        variable ``flags`` with the attributes ``long_name``, ``flag_names``, ``flag_descriptions`` and
        the statistics ``flag_<name>``.

        Parameters
        ----------
        output : {'prod', 'flags', None}
            - None: return nothing; the dataset of the masking object is updated.
            - ``'prod'``: return the full dataset updated with the variable ``flags``.
            - ``'flags'``: return the ``flags`` DataArray (its attributes include the flag statistics).

        Returns
        -------
        xarray.Dataset or xarray.DataArray or None
            See `output`.
        '''

        self.flag_names[:] = None
        self.flag_descriptions[:] = None
        self.flag_stats = {}

        # Other masks based on dask and tile approach
        steps = [self._nodata_step()]

        steps.append(self._s2cloudless_step() if self.sentinel2 else self._landsat_step())
        steps.append(self._water_step())

        if self.bcirrus:
            steps.append(self._cirrus_step())
        steps.append(self._high_swir_step())

        steps.append(self._surfwater_step())

        flags_lazy = self._flags_lazy(steps)
        bits = self._used_bits()
        counts_lazy = [((flags_lazy & np.uint32(1 << bit)) != 0).sum() for bit in bits]

        # one evaluation of the graph for the flags and for all the statistics
        flags, *counts = dask.compute(flags_lazy, *counts_lazy)

        self.flags = flags

        # OmniCloudMask on full raster
        self._get_omnicloudmask()

        self._store_stats(bits, counts, flags.size)

        # no copy: the dataset variable wraps the same uint32 array
        self.prod[FLAG_NAME] = (('y', 'x'), self.flags)
        self.prod[FLAG_NAME].attrs['long_name'] = 'flags computed from l1c image'
        self.prod[FLAG_NAME].attrs['flag_descriptions'] = self.flag_descriptions.astype(str)
        self.prod[FLAG_NAME].attrs['flag_names'] = self.flag_names.astype(str)

        if output == "prod":
            return self.prod
        elif output == "flags":
            flags = self.prod[FLAG_NAME]
            flags.attrs.update(self.flag_stats)
            return flags

    @staticmethod
    def _mask_block(block, tomask_value, tokeep_value, _type):
        '''Bad-pixel mask of a block of flags (numpy); see `create_mask` for the rules.'''
        if tomask_value and tokeep_value:
            bad = ((block & tomask_value) != 0) | ((block & tokeep_value) == 0)
        elif tokeep_value:
            bad = (block & tokeep_value) == 0
        elif tomask_value:
            bad = (block & tomask_value) != 0
        else:
            return np.zeros(block.shape, dtype=_type)
        return bad.astype(_type)

    @staticmethod
    def create_mask(flags,
                    tomask=(0, 2),
                    tokeep=(3,),
                    mask_name="mask",
                    _type=np.uint8
                    ):
        '''
        Create a binary mask from bit flags, with selection of the bits to mask or to keep.

        A pixel is bad (mask == 1) if **any** bit of `tomask` is set, **or** if **none** of the bits of
        `tokeep` is set (so with several `tokeep` bits a pixel is kept when at least one is set).
        Otherwise it is good (mask == 0). Empty `tomask` / `tokeep` disable the corresponding test.

        Memory: numpy-backed flags are processed by blocks of rows straight into the (1 byte/pixel)
        output; dask-backed flags give a lazy mask computed chunk by chunk (parallel).

        Parameters
        ----------
        flags : xarray.DataArray
            Bit flags (integer).
        tomask : sequence of int
            Bit numbers that make a pixel bad.
        tokeep : sequence of int
            Bit numbers of which at least one must be set for a pixel to be kept.
        mask_name : str
            Name of the output mask.
        _type : numpy dtype
            dtype of the mask (uint8 is recommended).

        Returns
        -------
        xarray.DataArray
            The mask (0 = good, 1 = bad).

        Examples
        --------
        >>> mask = Masking.create_mask(raster.flags,
        ...                            tomask=[0, 2, 11],
        ...                            tokeep=[3],
        ...                            mask_name="mask_from_flags")
        '''
        tomask_value = sum(1 << b for b in tomask)
        tokeep_value = sum(1 << b for b in tokeep)

        if flags.chunks is not None:  # dask-backed: stay lazy
            mask = xr.apply_ufunc(Masking._mask_block, flags,
                                  kwargs=dict(tomask_value=tomask_value, tokeep_value=tokeep_value, _type=_type),
                                  dask='parallelized', output_dtypes=[_type])
        else:
            values = flags.values
            out = np.empty(values.shape, dtype=_type)
            rows = max(1, (1 << 22) // max(values.shape[-1], 1))
            for y0 in range(0, values.shape[0], rows):
                out[y0:y0 + rows] = Masking._mask_block(values[y0:y0 + rows], tomask_value, tokeep_value, _type)
            mask = xr.DataArray(out, coords=flags.coords, dims=flags.dims)
        mask.attrs["long_name"] = "binary mask from flags"
        mask.attrs["description"] = "good pixels for mask == 0, bad pixels when mask == 1"
        mask.name = mask_name
        return mask
