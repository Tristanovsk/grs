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

=====  ==============================  ==========================================================================
bit    name                            meaning
=====  ==============================  ==========================================================================
0      ``nodata``                      NaN or crazy value (> 1e3) in the input
1-2    (unused)                        formerly s2cloudless cloud masks (GRS < 3.1), replaced by OmniCloudMask
3      ``water_swir_visible_index``    (visible - SWIR) / (visible + SWIR) above threshold
4      ``water_red_visible_index``     (visible - NIR) / (visible + NIR) above threshold (NDWI-like)
5      ``thin_cirrus``                 cirrus band above the thin cirrus threshold
6      ``opac_cirrus``                 cirrus band above the opaque cirrus threshold
7      ``high_swir``                   SWIR2 reflectance above threshold (bright cloud, glint, ...)
8      ``surfwater_land``              input variable ``surfwater`` == 0
9      ``surfwater_water``             input variable ``surfwater`` == 1
10     ``surfwater_cloud_and_shadow``  input variable ``surfwater`` >= 2
11     ``OCM_Thick_Cloud``             OmniCloudMask thick cloud (run at ``OCM_resolution``, 60 m by default)
12     ``OCM_Thin_Cloud``              OmniCloudMask thin cloud
13     ``OCM_Cloud_Shadow``            OmniCloudMask cloud shadow
14-20  ``l1_*``                        Landsat level-1 QA bits copied from ``flags_l1``  [Landsat only]
=====  ==============================  ==========================================================================

Clear sky
=========

The cloud and cloud shadow screening relies on OmniCloudMask (Sentinel-2 and Landsat): a pixel is clear
when its OmniCloudMask class is 0, i.e. none of the bits ``OCM_BITS`` (11, 12, 13) is set
(:meth:`Masking.clear_sky`). The processing chain only processes the clear water pixels (see
:py:meth:`grs.grs_process.Process.execute`).

How it is computed (dask)
=========================

* The flags are built as **one lazy dask graph**: the band cube is cut into ``tile_size`` x
  ``tile_size`` tiles (all wavelengths in one chunk) and every tile goes through all the masks in
  a single task. Tiles are processed in parallel by the active dask scheduler (threads by default,
  or a ``dask.distributed`` cluster if a client exists); there is no Python tile loop.
* Only the bands a mask needs are read (band selection is fused into the read when the dataset
  chunks match ``tile_size``).
* Neighbourhood operations (steps with a ``halo``) use overlapping tiles (``dask.array.overlap``)
  with a halo larger than their radius, so the result is identical to a run on the full image.
* :meth:`Masking.process` evaluates the tile flags **and** their statistics in one ``dask.compute``
  call, so every tile is processed exactly once; OmniCloudMask then runs once on a downscaled image
  (``OCM_resolution``). The result (``uint32``, 4 bytes per pixel) is kept
  in memory and stored in the dataset: cloud detection is never repeated downstream.
* Peak memory is about ``n_workers x tile_size**2 x (n_bands x 4 bytes x ~3)``. Lower
  ``tile_size`` or the number of workers (``dask.config.set(num_workers=...)``) if memory is tight.
'''

import logging
from collections import namedtuple
from importlib.metadata import version, PackageNotFoundError

import dask
import dask.array as da
from dask.array.overlap import overlap, trim_internal
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd

import logging

from s2cloudless import S2PixelCloudDetector

FLAG_NAME = 'flags'
# flag bits of the OmniCloudMask classes 1 (thick cloud), 2 (thin cloud) and 3 (cloud shadow)
OCM_BITS = (11, 12, 13)


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
    OCM_Rband, OCM_Gband, OCM_NIRband : float
        Red, green and NIR bands (nm, nearest available) given to OmniCloudMask.
    OCM_resolution : float or None
        Resolution (m) at which OmniCloudMask runs (default 60; e.g. 100 for speed); the mask is brought
        back to the native grid. None, or a value not above the native pixel size, runs at native resolution.
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

        # settings OmniCloudMask
        self.OCM_Rband = 665
        self.OCM_Gband = 560
        self.OCM_NIRband = 865
        self.OCM_batch_size = 1
        self.OCM_inference_dtype = "fp32"
        self.OCM_inference_device = 'cpu'
        self.OCM_mosaic_device = 'cpu'
        # resolution (m) at which OmniCloudMask runs (block-averaged bands); the mask is then brought
        # back to the native grid. None or <= native pixel size: run at native resolution
        self.OCM_resolution = 60

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
        as the largest halo of the steps (none of the current steps needs one).

    Attributes
    ----------
    flags : numpy.ndarray or None
        ``uint32`` flags (y, x), available after `process` or any single-mask method.
    flag_names, flag_descriptions : numpy.ndarray of object
        Name and description of every bit (None for unused bits).
    flag_stats : dict
        Fraction of flagged pixels per flag, keys ``'flag_<name>'``.
    sentinel2 : bool
        True for Sentinel-2, False for Landsat (the level-1 QA bits are copied).
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

    def nodata_mask(self,
                    bitmask=0,
                    name='nodata',
                    description='nodata in input image'):
        '''
        set flag nodata: condition either nan or crazy numerical values

        :return:
        '''

        self.flags = (np.isnan(self.prod.bands.isel(wl=1)) | (self.prod.bands.isel(wl=1) > 1e3)).astype(np.int32)
        # set naming and description as attributes
        self.flag_descriptions[bitmask] = description
        self.flag_names[bitmask] = name

    def cloud_mask(self,
                   bitmasks=[1, 2],
                   names=['cloud_p06', 'cloud_p08'],
                   descriptions=['low confidence cloud mask from s2cloudless with settings proba.',
                                 'high confidence cloud mask from s2cloudless with settings proba.']
                   ):
        '''
        Apply s2cloudless masking with two levels of confidence

        :param bitmasks: bit number on which the flag is coded
        :param names: name of the flag
        :param descriptions: description of the flag
        :return:
        '''

        logging.info('cloud masking with s2cloudless')

        # apply cloud mask
        bands = np.array([self.prod.bands.values.transpose((1, 2, 0))])
        # ---------------------
        # first low confidence
        # ---------------------
        bitmask = bitmasks[0]

        cloud_detector = S2PixelCloudDetector(threshold=self.low_confid_cloud_proba_thresh, average_over=1,
                                              dilation_size=self.low_confid_cloud_dilation, all_bands=True)
        probability_maps = cloud_detector.get_cloud_probability_maps(bands)
        cloud_mask = cloud_detector.get_mask_from_prob(probability_maps)[0]
        self.flags = self.flags + ((cloud_mask == 1) << bitmask)

        # add name and description
        self.flag_descriptions[bitmask] = descriptions[0] + ' threshold={:.2f}, dilation_size={:d}'.format(
            self.low_confid_cloud_proba_thresh, self.low_confid_cloud_dilation)
        self.flag_names[bitmask] = names[0]

        # ---------------------
        # second high confidence
        # ---------------------
        bitmask = bitmasks[1]

        cloud_detector = S2PixelCloudDetector(threshold=self.high_confid_cloud_proba_thresh, average_over=1,
                                              dilation_size=self.high_confid_cloud_dilation, all_bands=True)
        cloud_mask = cloud_detector.get_mask_from_prob(probability_maps)[0]
        self.flags = self.flags + ((cloud_mask == 1) << bitmask)

        # add name and description
        self.flag_descriptions[bitmask] = descriptions[1] + ' threshold={:.2f}, dilation_size={:d}'.format(
            self.high_confid_cloud_proba_thresh, self.high_confid_cloud_dilation)
        self.flag_names[bitmask] = names[1]

        del cloud_detector, cloud_mask, probability_maps

    def water_mask(self,
                   bitmasks=[3, 4],
                   names=['water_swir_visible_index', 'water_red_visible_index'],
                   descriptions=['water mask from normalized index from swir and visible band',
                                 'water mask from normalized index from nir and visible band, warning could fail for turbid waters'],
                   ):
        '''
        apply water/land mask
        compute mask from NDWI from visible and NIR and visible and SWIR

        :param bitmasks: bit numbers on which the flags are coded
        :param names: name of the flags
        :param descriptions: description of the flags
        :return:
        '''

        logging.info('water masking')
        visible = self.prod.bands.sel(wl=self.bvis, method='nearest')
        nir = self.prod.bands.sel(wl=self.bnir, method='nearest')
        swir = self.prod.bands.sel(wl=self.bswir, method='nearest')

        ndwi = (visible - nir) / (visible + nir)
        ndwi_swir = (visible - swir) / (visible + swir)

        # set flags raster
        self.flags = self.flags + (
                ((ndwi_swir.values > self.vis_swir_index_threshold) << bitmasks[0]) +
                ((ndwi.values > self.ndwi_threshold) << bitmasks[1]))

        # set naming and description as attributes
        self.flag_descriptions[bitmasks[0]] = (descriptions[0] +
                                               ', bands centered on {:.1f} and {:.1f} nm'.format(self.bvis,
                                                                                                 self.bswir) +
                                               ', threshold={:.2f}'.format(self.vis_swir_index_threshold))
        self.flag_names[bitmasks[0]] = names[0]

        self.flag_descriptions[bitmasks[1]] = (descriptions[1] +
                                               ', bands centered on {:.1f} and {:.1f} nm'.format(self.bvis, self.bnir) +
                                               ', threshold={:.2f}'.format(self.ndwi_threshold))
        self.flag_names[bitmasks[1]] = names[1]

    def cirrus_mask(self,
                    bitmasks=[5, 6],
                    names=['thin_cirrus', 'opac_cirrus'],
                    descriptions=['thin cirrus mask from cirrus band',
                                  'opac cirrus mask from cirrus band']
                    ):
        '''
        Compute cirrus mask for thin and opac high clouds from spectral band around 1375 nm, if exists

        :param bitmasks: bit numbers on which the flags are coded
        :param names: name of the flags
        :param descriptions: description of the flags
        :return:
        '''

        logging.info('cirrus masking')
        cirrus = self.prod.bands.sel(wl=self.bcirrus, method='nearest')
        thresholds = [self.thin_cirrus_threshold, self.opac_cirrus_threshold]
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
            # set naming and description as attributes
            self.flag_descriptions[bitmasks[ii]] = descriptions[ii]
            self.flag_names[bitmasks[ii]] = names[ii]

    def duplicate_landsat_flags(self,
                                bitmasks=[11, 12, 13, 14, 15, 16, 17],
                                landsat_bitmasks=[1, 2, 3, 4, 5, 6, 7],
                                names=['l1_dilated_cloud', 'l1_cirrus', 'l1_cloud',
                                       'l1_cloud_shadow', 'l1_snow', 'l1_clear', 'l1_water'],
                                descriptions=['landsat dilated cloud mask',
                                              'landsat cirrus mask',
                                              'landsat cloud mask',
                                              'landsat cloud shadow mask',
                                              'landsat snow mask',
                                              'landsat clear mask',
                                              'landsat water mask'],
                                flags_landsat_band='flags_l1'
                                ):
        '''
        Save and duplicate landsat main masks from level 1 flags

        :param bitmasks: bit numbers on which the flags are coded
        :param landsat_bitmasks: bit numbers of the landsat flags
        :param names: name of the flags
        :param descriptions: description of the flags
        :param flags_landsat_band: name of the landsat flags band
        :return:
        '''

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
        Generate the flags raster and attributes

        :param output:
            - if None returns nothing but the raster is updated within the masking object
            - if "prod" returns the full raster updated with the flags variable and attributes
            - if "flags" returns xarray DataArray of the flags
              plus the dictionary of flags statistics
        :return: see :param output
        '''
        # apply the masking processors
        self.nodata_mask()
        if self.sentinel2:
            self.cloud_mask()
        else:
            self.duplicate_landsat_flags()
        self.water_mask()
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

        # OmniCloudMask on a downscaled image, mask brought back to the native grid
        counts += self._get_omnicloudmask(bitmasks=OCM_BITS)
        bits = bits + list(OCM_BITS)

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
    def clear_sky(flags):
        '''
        Clear-sky pixels: OmniCloudMask class 0 (no thick cloud, thin cloud or cloud shadow), i.e. none of
        the bits ``OCM_BITS`` set.

        :param flags: ``flags`` DataArray (numpy- or dask-backed)
        :return: boolean DataArray, True for the clear pixels
        '''
        return (flags & np.uint32(sum(1 << b for b in OCM_BITS))) == 0

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
