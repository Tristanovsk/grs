'''
Module to set output product and format (L2A): NetCDF and Zarr export.

An :class:`L2aProduct` gathers the L2A variables (remote-sensing reflectance, aerosol, geometry, flags,
DEM, ...) and the coarse-resolution ancillary data (CAMS fields, gaseous transmittance), and writes them
with a fixed, lossy-but-controlled packing (``int16`` with scale/offset, see `_base_encoding`).

Choosing the format from the file name
======================================

``export(output_path)`` picks the format from the extension of `output_path`:

=================  ========  =========================================================
extension          format    files
=================  ========  =========================================================
``.nc``            NetCDF    ``<output_path>`` (main) and ``<stem>_anc.nc`` (ancillary)
``.zarr``          Zarr      ``<output_path>`` (store, with pyramids)
none               NetCDF    legacy: ``<output_path>/<basename>.nc`` and ``_anc.nc``
                             (or the format given by ``fmt``)
=================  ========  =========================================================

``fmt='netcdf'`` / ``fmt='zarr'`` forces the format; with an extension that says otherwise a
``ValueError`` is raised. Without extension and ``fmt``, ``output_path`` is a directory whose name is the
base name of the files (previous behaviour). With ``fmt='zarr'`` and no extension the store is
``<output_path>/<basename>.zarr``.

Output layouts
==============

NetCDF (unchanged layout)::

    <dir>/<stem>.nc        main product
    <dir>/<stem>_anc.nc    ancillary data (coarse grid, dimensions xc, yc)

Zarr (with pyramids, see below)::

    <name>.zarr/                  root group: product attributes and ``multiscales``, no data
    <name>.zarr/0/                main product at full resolution (level 0)
    <name>.zarr/1/                same variables, 2x coarser
    <name>.zarr/2/                4x coarser
    ...
    <name>.zarr/ancillary/        ancillary data (coarse grid, dimensions xc, yc)

Read back with ``xr.open_zarr(store, group='0')`` (full resolution), ``xr.open_zarr(store, group='2')`` and
``xr.open_zarr(store, group='ancillary')``. The packed variables are unpacked automatically
(``scale_factor`` / ``add_offset`` / ``_FillValue``). The root attribute ``multiscales`` lists the levels
(group, scale factor, shape) for viewers.

Pyramids
========

Each level is built from the previous one (factor 2 per level, so level ``k`` is ``2**k`` coarser):

* continuous variables (``Rrs``, ``aot550``, angles, ``dem``, ...): mean of the 2 x 2 block, ignoring NaN;
* integer variables (``flags``, ``mask``): value of the top-left pixel of the block (nearest). Averaging or
  OR-ing bit flags would create meaningless values; the overviews are for display, use the full resolution
  for any statistics on the flags;
* the coordinates are the block centres and the georeferencing (``GeoTransform``) is updated per level.

The number of levels is chosen so that the coarsest one is at most ``pyramid_min_size`` pixels (512) on
its longest side (or set with ``pyramid_levels``). The overviews add about 1/3 of the full-resolution size
(1/4 + 1/16 + ...). ``pyramids=False`` writes the full resolution only.

Cloud-ready Zarr
================

The store is laid out for reading a window of the image without downloading the rest, and for being
written straight to object storage:

* **Sharding** (zarr-python >= 3, Zarr format 3, and an xarray with the ``shards`` encoding): small inner
  chunks (default 256 x 256 pixels, all wavelengths) are packed into large shards (default 1024 x 1024).
  A reader fetches one inner chunk with a byte-range request, while the store only holds a few large
  objects (4x4 = 16 chunks per object) instead of thousands of small ones: fewer requests to list / write,
  no per-object overhead. The shard shape is the dask chunk and a multiple of the inner chunk; a level
  that is a single block (small overviews) is stored without sharding. If a requirement is missing the
  export warns and writes plain 512 x 512 chunks (Zarr format 3 with zarr-python 3, format 2 otherwise).
* **Compression**: Blosc/zstd with byte shuffle (``compression`` / ``clevel``), well suited to the packed
  ``int16`` data and the bit flags.
* **Consolidated metadata**: the metadata of the whole hierarchy (full resolution, overviews, ancillary)
  is read in one request when opening the store.
* **Overviews** (see Pyramids) for fast low-resolution display.
* **Remote output**: ``output_path`` can be an fsspec URL (``s3://bucket/scene.zarr``, ``gs://``,
  ``az://``...; credentials through ``storage_options``). The chunks are written directly to the bucket as
  they are computed, nothing is staged on the local disk. Needs the matching fsspec driver (``s3fs``,
  ``gcsfs``, ``adlfs``). NetCDF is local only.

Dask chunks are the shards, so every shard is written by exactly one task (concurrent writes into the same
shard would corrupt it).

Zarr export and dask
====================

The product is normally dask-backed (lazy). The Zarr export is built to evaluate it **once**:

* the main and ancillary datasets are rechunked to regular chunks (default ``y = x = 1024`` and all
  wavelengths in one chunk, the same layout ``Product`` reads) so Zarr chunks and dask chunks coincide;
* each overview is a lazy coarsening of the previous level, so it reuses the chunks of the full resolution
  graph instead of recomputing them;
* all the writes (full resolution, overviews, ancillary) are prepared with ``compute=False`` and run
  together in one ``dask.compute`` call, so intermediate results shared by the graphs (e.g. the gaseous
  transmittance, the full-resolution chunks) are computed only once and chunks are written in parallel as
  soon as they are ready;
* the packing is fixed (no data-dependent range), so no extra pass over the data is needed;
* the metadata of the whole hierarchy is consolidated at the end (one small file read to open it,
  useful on object storage).

Neither export modifies the product: the SNAP-compliant band explosion, the attribute clean-up and
the rechunking are done on a shallow copy, so ``export_to_netcdf`` and ``export_to_zarr`` can be called
one after the other on the same object.
'''
import os
import shutil
import numpy as np
import xarray as xr
import logging


import importlib_resources
import yaml

configfile = importlib_resources.files(__package__) / 'config.yml'
with open(configfile, 'r') as file:
    config = yaml.safe_load(file)

NETCDF_ENGINE = 'netcdf4' #config['processor']['netcdf_engine']

# attributes that would conflict with the encoding when written
_CONFLICTING_ATTRS = ('grid_mapping', 'scale_factor', '_FillValue', 'add_offset')

# file extension -> format
_NETCDF_EXTENSIONS = ('.nc', '.nc4', '.cdf')
_ZARR_EXTENSIONS = ('.zarr',)
_FORMAT_EXTENSION = {'netcdf': '.nc', 'zarr': '.zarr'}

# pyramids
PYRAMID_FACTOR = 2
PYRAMID_MIN_SIZE = 512
LEVEL_GROUP = '{:d}'
ANCILLARY_GROUP = 'ancillary'

# Zarr layout (pixels; -1 = whole dimension): write blocks (= shards), stored chunks
_DEFAULT_SHARDS = {'y': 1024, 'x': 1024, 'wl': -1}
_DEFAULT_CHUNKS_SHARDED = {'y': 256, 'x': 256, 'wl': -1}
_DEFAULT_CHUNKS = {'y': 512, 'x': 512, 'wl': -1}


def _base_encoding():
    '''
    Fixed packing of the L2A variables, independent of the output format
    (compression / chunking options are added by the exporters).

    ============  =======  =============  ==========  ============
    variable      dtype    scale_factor   add_offset  _FillValue
    ============  =======  =============  ==========  ============
    aot550        int16    0.001          -           -9999
    BRDFg         int16    1e-5           0.3         -32768
    dem           int16    1              500         -32768
    vza           int16    0.01           -           -32768
    raa           int16    0.01           180         -32768
    sza           int16    0.01           30          -32768
    Rrs           int16    1e-5           0.3         -32768
    flags, mask   (unchanged dtype, no packing)
    ============  =======  =============  ==========  ============

    :return: new dictionary ``{variable: encoding}``
    '''
    gm = 'spatial_ref'
    return {
        'aot550': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, 'grid_mapping': gm},
        'BRDFg': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768,
                  'grid_mapping': gm},
        'dem': {'dtype': 'int16', 'scale_factor': 1, 'add_offset': 500, '_FillValue': -32768,
                'grid_mapping': gm},
        'vza': {'dtype': 'int16', 'scale_factor': 0.01, '_FillValue': -32768, 'grid_mapping': gm},
        'raa': {'dtype': 'int16', 'scale_factor': 0.01, 'add_offset': 180, '_FillValue': -32768,
                'grid_mapping': gm},
        'sza': {'dtype': 'int16', 'scale_factor': 0.01, 'add_offset': 30, '_FillValue': -32768,
                'grid_mapping': gm},
        'flags': {'grid_mapping': gm},
        'mask': {'grid_mapping': gm},
    }


def _rrs_encoding():
    '''Packing of the remote-sensing reflectance (and its per-band variants): int16, 1e-5, offset 0.3.'''
    return {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3,
            '_FillValue': -32768, 'grid_mapping': 'spatial_ref'}


def _coord_encoding():
    '''
    Encoding of the non-dimension coordinates (scalars): ``time`` as CF integer seconds,
    ``spatial_ref`` (CRS container, value meaningless) as a small integer.

    :return: new dictionary ``{coordinate: encoding}``
    '''
    return {
        'time': {'units': 'seconds since 1970-01-01T00:00:00', 'calendar': 'proleptic_gregorian',
                 'dtype': 'int64'},
        'spatial_ref': {'dtype': 'int32'},
    }


def _package_version(name):
    '''Installed version of the distribution `name` (from its metadata, without importing it); 'unknown' if absent.'''
    from importlib.metadata import version, PackageNotFoundError
    try:
        return version(name)
    except PackageNotFoundError:
        return 'unknown'


def _reset_encoding(ds):
    '''
    Drop, in place, the encodings inherited from the input files (chunks, codecs, time units...) of all the
    variables of `ds`, coordinates included, so they cannot conflict with the Zarr layout / version.
    `ds` must be a shallow copy.
    '''
    for var in ds.variables:
        ds[var].encoding = {}


def _load_coords(ds):
    '''
    Dataset with its non-dimension coordinates (``time``, ``spatial_ref``: scalars) held in memory.

    Lazy (dask) coordinates, e.g. read from the L1C file or chunked by ``Dataset.chunk``, would otherwise be
    written by the delayed ``to_zarr`` of every group inside the final ``dask.compute``; in memory they are
    written with the metadata.
    '''
    lazy = {c: ds[c].variable.compute() for c in ds.coords if c not in ds.dims and ds[c].chunks is not None}
    return ds.assign_coords(lazy) if lazy else ds


def _format_from_name(path):
    '''
    Output format deduced from the extension of `path`.

    :param path: output path
    :return: ``'netcdf'`` (``.nc``, ``.nc4``, ``.cdf``), ``'zarr'`` (``.zarr``) or None (other / no extension)
    '''
    ext = os.path.splitext(str(path).rstrip('/\\'))[1].lower()  # no normpath: it would break 's3://...'
    if ext in _NETCDF_EXTENSIONS:
        return 'netcdf'
    if ext in _ZARR_EXTENSIONS:
        return 'zarr'
    return None


def _is_remote(path):
    '''True if `path` is an fsspec URL (``s3://``, ``gs://``, ``az://``, ``https://``...), False for a local path.'''
    path = str(path)
    return '://' in path and not path.startswith('file://')


def _zarr_major():
    '''Major version of zarr-python, None if zarr is not installed.'''
    try:
        import zarr
        return int(zarr.__version__.split('.')[0])
    except Exception:
        return None


def _xarray_supports_shards():
    '''True if the installed xarray accepts the ``shards`` key in the Zarr encoding (feature test, not a version test).'''
    try:
        import inspect
        import xarray.backends.zarr as xarray_zarr
        return 'shards' in inspect.getsource(xarray_zarr)
    except Exception:
        return False


def _compression_encoding(compression, clevel, zarr_major):
    '''
    Compression options for the encoding of the Zarr variables: Blosc with `compression` codec
    (``'zstd'``, ``'lz4'``, ``'zlib'``...) and byte shuffle. The option name depends on the zarr-python
    version (``compressors`` for 3, ``compressor`` for 2).

    :return: dict to merge in the encoding of each variable ({} if `compression` is None or unavailable)
    '''
    if compression is None or zarr_major is None:
        return {}
    try:
        if zarr_major >= 3:
            from zarr.codecs import BloscCodec
            return {'compressors': (BloscCodec(cname=compression, clevel=clevel, shuffle='shuffle'),)}
        import numcodecs
        return {'compressor': numcodecs.Blosc(cname=compression, clevel=clevel, shuffle=numcodecs.Blosc.SHUFFLE)}
    except Exception as e:
        logging.warning(f'compression {compression!r} not available, using the Zarr default: {e}')
        return {}


def _var_chunking(dims, dask_chunks, chunks):
    '''
    Stored chunk shape and write-block (shard) shape of a variable.

    The block is the dask chunk (regular chunks: first one). The stored chunk is the requested one, clipped
    to the block; it has to divide the block, otherwise the whole block is used (small overviews).

    :param dims: dimension names of the variable
    :param dask_chunks: ``variable.chunks`` (tuple of tuples)
    :param chunks: requested stored chunk size per dimension (-1 / missing: whole block)
    :return: ``(chunk_shape, block_shape)``, tuples in the order of `dims`
    '''
    inner, block = [], []
    for dim, dc in zip(dims, dask_chunks):
        blk = int(dc[0])
        want = chunks.get(dim, -1)
        c = blk if (want is None or want < 0) else min(int(want), blk)
        if blk % c:
            c = blk
        inner.append(c)
        block.append(blk)
    return tuple(inner), tuple(block)


def _clear_store(store, storage_options=None):
    '''Delete an existing local or remote Zarr store (a store left from a previous run).'''
    if _is_remote(store):
        import fsspec
        fs, path = fsspec.core.url_to_fs(store, **(storage_options or {}))
        if fs.exists(path):
            fs.rm(path, recursive=True)
    elif os.path.exists(store):
        shutil.rmtree(store)


def _consolidate(store, storage_options=None, zarr_major=None, groups=()):
    '''
    Consolidate the metadata of the local or remote store: the root group (which covers the whole
    hierarchy) and each group of `groups` on its own, so that a group opened directly
    (``xr.open_zarr(store, group='2')``) also needs a single metadata request.

    :param store: path or URL of the store
    :param storage_options: fsspec options for a URL
    :param zarr_major: major version of zarr-python
    :param groups: paths of the groups to consolidate individually (e.g. ``['0', '1', 'ancillary']``)
    '''
    import zarr
    if _is_remote(store):
        if (zarr_major or 0) >= 3:
            # zarr 3 needs an async filesystem (s3fs, gcsfs, adlfs...)
            target = zarr.storage.FsspecStore.from_url(store, storage_options=storage_options or {})
        else:
            import fsspec
            target = fsspec.get_mapper(store, **(storage_options or {}))
    else:
        target = store
    zarr.consolidate_metadata(target)
    for group in groups:
        zarr.consolidate_metadata(target, path=group)


def _pyramid_levels(ny, nx, min_size=PYRAMID_MIN_SIZE, max_levels=None):
    '''
    Number of overview levels (not counting the full resolution): the image is halved (rounding up) until
    its longest side is at most `min_size`.

    >>> _pyramid_levels(5490, 5490)
    4
    >>> _pyramid_levels(400, 300)
    0

    :param ny, nx: size of the full resolution image
    :param min_size: size under which the coarsest level has to be (pixels, longest side)
    :param max_levels: upper bound of the number of levels (None: no bound)
    :return: number of overviews (>= 0)
    '''
    levels = 0
    while max(ny, nx) > min_size and min(ny, nx) > 1:
        ny, nx = -(-ny // PYRAMID_FACTOR), -(-nx // PYRAMID_FACTOR)  # ceil
        levels += 1
        if max_levels is not None and levels >= max_levels:
            break
    return levels


def _coarsen(ds, factor=PYRAMID_FACTOR):
    '''
    Lazy reduction of the spatial resolution of a dataset by `factor` (dims ``y`` and ``x``).

    Continuous (float) variables: mean of the blocks ignoring NaN. Integer / boolean variables (flags,
    masks): value of the first pixel of the block (nearest). Incomplete blocks at the right / bottom edges
    are kept (padded), so the size is ``ceil(n / factor)``. The coordinates are the block centres, the
    georeferencing is updated.

    :param ds: xarray.Dataset with dimensions ``y`` and ``x``
    :param factor: integer reduction factor
    :return: coarser xarray.Dataset (dask-backed if `ds` is)
    '''
    window = {'y': factor, 'x': factor}
    is_spatial = [v for v in ds.data_vars if 'y' in ds[v].dims and 'x' in ds[v].dims]
    discrete = [v for v in is_spatial if ds[v].dtype.kind in 'biu']
    continuous = [v for v in ds.data_vars if v not in discrete]

    coarse = ds[continuous].coarsen(window, boundary='pad', coord_func='mean').mean(keep_attrs=True)
    if discrete:
        nearest = ds[discrete].isel(y=slice(0, None, factor), x=slice(0, None, factor))
        nearest = nearest.assign_coords(y=coarse['y'], x=coarse['x'])
        coarse = coarse.assign({v: nearest[v] for v in discrete})
    coarse.attrs = dict(ds.attrs)

    # geotransform of the coarse grid, recomputed from the coordinates (first pixel = block centre)
    return _update_geotransform(coarse)


def _update_geotransform(ds, dims=('y', 'x')):
    '''
    Rewrite the ``GeoTransform`` of the ``spatial_ref`` variable from the coordinates of `dims`.

    Needed after resampling (overviews, coarse ancillary grid) so that GIS readers see the right pixel size
    and origin. Failures (e.g. no CRS or non-georeferenced data) only log a warning: georeferencing is
    informative only and never blocks the export.

    :param ds: xarray.Dataset
    :param dims: names of the (y, x) spatial dimensions
    :return: dataset with the updated ``GeoTransform`` (unchanged on failure)
    '''
    try:
        ds = ds.rio.set_spatial_dims(x_dim=dims[1], y_dim=dims[0])
        ds = ds.rio.write_transform(ds.rio.transform(recalc=True))
    except Exception as e:
        logging.warning(f'could not update the geotransform: {e}')
    return ds


class L2aProduct():
    '''
    Create and handle the L2A product object.

    Parameters
    ----------
    prod : object
        Native product (``Product``-like): ``prod.raster`` is the dataset of the input image (with
        ``vza``, ``sza``, ``raa``, ``flags``, ``wl_true``, attributes), and ``prod`` carries the processing
        parameters ``sunglint_threshold``, ``ndwi_threshold``, ``vis_swir_index_threshold``,
        ``hcld_threshold``, ``dirdata``, ``abs_gas_file``, ``water_vapor_transmittance_file`` (stored as
        string attributes of the output) and optionally ``bcirrus``/``cirrus``, ``bwv``/``wv``.
    l2_prod : xarray.Dataset
        L2A variables computed by the processor (e.g. ``Rrs``, ``aot550``, ``BRDFg``), on the image grid.
    cams : object
        CAMS ancillary data, with ``cams.raster`` (dataset with dimensions ``x``, ``y``).
    gas_trans : object
        ``GaseousTransmittance`` with the coarse transmittance ``Tg_tot_coarse`` (dimensions ``x``, ``y``, ``wl``).
    dem : xarray.DataArray, optional
        Digital elevation model on the image grid; a NaN-filled ``dem`` variable is written if None.

    Attributes
    ----------
    l2_prod : xarray.Dataset
        Main product after `construct_l2a` (merged with ``flags`` and ``dem``, georeferenced).
    ancillary : xarray.Dataset
        Coarse-resolution ancillary data (dimensions ``xc``, ``yc``; ``wl`` renamed ``wl_all``).
    complevel : int
        zlib compression level used for NetCDF (default 5).
    wl_cirrus_band, wl_wv_band : int
        Wavelengths (nm) of the cirrus (1375) and water vapour (945) bands.
    '''

    def __init__(self, prod, l2_prod, cams, gas_trans, dem=None):
        self.prod = prod
        self.l2_prod = l2_prod
        self.cams = cams
        self.gas_trans = gas_trans
        self.dem = dem
        self.ancillary = None
        self.wl_cirrus_band = 1375
        self.wl_wv_band = 945
        self.complevel = 5
        self.construct_l2a()

    def construct_l2a(self):
        '''
        Assemble the main and the ancillary datasets with their metadata and georeferencing.

        Main product (``self.l2_prod``): attributes of the native raster plus the processing parameters
        (as strings) and ``GRSdriver_version``, the geometry angles (``vza``, ``raa`` averaged over
        wavelengths if band dependent; ``sza``), ``flags`` and ``dem``; spatial dimensions ``x``/``y``, CRS written.
        Ancillary (``self.ancillary``): coarse transmittance and CAMS fields on dimensions ``xc``/``yc``.

        Notes
        -----
        * The attributes are copied, so the native raster's attributes are no longer modified.
        * ``cirrus_band`` and ``wv_band`` are built from ``prod.cirrus`` / ``prod.wv`` but are not merged
          into the output (as in the previous version); add them to ``l2_prod`` if they should be exported.
        * The rioxarray calls no longer use ``inplace=True`` (removed in recent rioxarray versions); their
          results are assigned instead.
        '''
        logging.info('construct l2a')
        native_raster = self.prod.raster
        # first get attributes from native product (copy: do not modify the native raster)
        self.l2_prod.attrs = dict(native_raster.attrs)

        keys = ['sunglint_threshold', 'ndwi_threshold', 'vis_swir_index_threshold', 'hcld_threshold',
                'dirdata', 'abs_gas_file', 'water_vapor_transmittance_file']

        for key in keys:
            self.l2_prod.attrs[key] = str(self.prod.__dict__[key])
        # version of the L1C reader
        self.l2_prod.attrs['GRSdriver_version'] = _package_version('GRSdriver')

        # rename dim for coarse resolution data
        cams_raster = self.cams.raster.rename({"x": "xc", "y": "yc"})
        transmittance_raster = self.gas_trans.Tg_tot_coarse.rename({"x": "xc", "y": "yc", 'wl': 'wl_all'})

        # save geometrie angles, take mean values if angle depends on spectral bands
        def angle_extraction(angle):
            if 'wl' in angle.coords:
                return angle.mean('wl')
            else:
                return angle

        self.l2_prod['vza'] = angle_extraction(native_raster.vza)
        self.l2_prod['sza'] = native_raster.sza
        self.l2_prod['raa'] = angle_extraction(native_raster.raa)

        # add cirrus and water vapor bands
        if self.prod.bcirrus:
            cirrus = self.prod.cirrus
            cirrus.name = 'cirrus_band'
            cirrus.attrs = {
                'description': 'TOA reflectance in cirrus band from L1C image, might be used to filter out high clouds',
            }

        if self.prod.bwv:
            wvband = self.prod.wv
            wvband.name = 'wv_band'
            wvband.attrs = {
                'description': 'TOA reflectance in water vapor band from L1C image',
            }

        # ------------------
        # adding digital elevation model to output
        dem = self.dem
        # add empty DEM in raster
        if dem is None:
            dem = xr.full_like(self.l2_prod['sza'], np.nan, dtype=np.float32)
            dem.name = 'dem'

        # final merge
        self.l2_prod = xr.merge([self.l2_prod, native_raster.flags, dem], compat='override')
        self.l2_prod = self.l2_prod.rio.set_spatial_dims(x_dim='x', y_dim='y')
        self.l2_prod = self.l2_prod.rio.write_coordinate_system()
        self.l2_prod = self.l2_prod.rio.write_crs()
        self.l2_prod = _update_geotransform(self.l2_prod)

        self.ancillary = xr.merge([transmittance_raster, cams_raster], compat='override')
        self.ancillary = self.ancillary.rio.set_spatial_dims(x_dim='xc', y_dim='yc')
        self.ancillary = self.ancillary.rio.write_coordinate_system()
        self.ancillary = self.ancillary.rio.write_crs()
        # coarse grid: the GeoTransform inherited from the inputs does not match its pixel size
        self.ancillary = _update_geotransform(self.ancillary, dims=('yc', 'xc'))

    # ------------------------------------------------------------------
    # common preparation of the export
    # ------------------------------------------------------------------
    @staticmethod
    def _scrub_attrs(ds):
        '''Remove, in place, the variable attributes that conflict with the encoding (``grid_mapping``,
        ``scale_factor``, ``_FillValue``, ``add_offset``). `ds` must be a shallow copy.'''
        for var in ds.keys():
            for param in _CONFLICTING_ATTRS:
                ds[var].attrs.pop(param, None)

    def _prepare_export(self, snap_compliant=False):
        '''
        Dataset and encoding to write, built on a shallow copy of the product (``self`` is not modified).

        :param snap_compliant: if True, ``Rrs`` is exploded into one variable per band (``Rrs_<wl>``) with
            the attributes expected by SNAP, and the dataset is flagged with ``metadata_profile = 'beam'``
            (``'datacube'`` otherwise)
        :return: (dataset, encoding) where the encoding only contains variables present in the dataset
            and no format-specific option (compression, chunks)
        '''
        ds = self.l2_prod.copy()
        encoding = _base_encoding()

        if snap_compliant:
            # explode Rrs array to get one variable per band to be SNAP compliant
            for var in ['Rrs']:
                img_snap = ds[var].to_dataset("wl")
                suff = ''
                if var == 'Rrs_g':
                    suff = 'with sunglint '
                for ii, band in enumerate(img_snap.keys()):
                    band_name = var + '_{:d}'.format(band)
                    wavelength = self.prod.raster.wl_true.values[ii]
                    img_snap[band].attrs = {
                        'long_name': 'Remote sensing reflectance' + suff + ' at {:d} nm'.format(band),
                        'Unit': 'sr-1',
                        'units': 'sr-1',
                        'radiation_wavelength': wavelength,
                        'radiation_wavelength_unit': 'nm',
                        'wavelength': wavelength,
                        'valid_pixel_expression': ''}
                    encoding[band_name] = _rrs_encoding()

                    img_snap = img_snap.rename({band: band_name})
                ds = xr.merge([ds.drop_vars(var), img_snap], compat='override')
            ds.attrs['auto_grouping'] = 'Rrs:Rrs_g'
            ds.attrs['metadata_profile'] = 'beam'
        else:
            encoding['Rrs'] = _rrs_encoding()
            ds.attrs['metadata_profile'] = 'datacube'

        # fix for conflicts with attrs and encoding, needs to remove 'grid_mapping' from input attrs
        self._scrub_attrs(ds)

        encoding = {k: v for k, v in encoding.items() if k in ds.variables}
        return ds, encoding

    @staticmethod
    def _output_names(output_path, ext=None):
        '''
        Resolve the output location and create its directory.

        * `output_path` has the extension `ext` (e.g. ``/data/scene.zarr``): it is the output file itself;
          the directory is its parent and the base name its stem (``scene``).
        * otherwise (legacy): `output_path` is a directory whose name is the base name; the file is
          ``<output_path>/<basename><ext>``.

        :param output_path: output path (robust to a trailing slash)
        :param ext: extension of the format (``'.nc'``, ``'.zarr'``), or None
        :return: ``(directory, basename, main_path)`` where `main_path` is the main file / store
            (equal to `directory`/`basename` without extension if `ext` is None)
        '''
        output_path = os.path.normpath(output_path)
        stem, current_ext = os.path.splitext(os.path.basename(output_path))
        if ext is not None and current_ext.lower() in (_NETCDF_EXTENSIONS + _ZARR_EXTENSIONS):
            directory = os.path.dirname(output_path) or '.'
            basename = stem
            main = output_path
        else:
            directory = output_path
            basename = os.path.basename(output_path)
            main = os.path.join(directory, basename + (ext or ''))
        os.makedirs(directory, exist_ok=True)
        return directory, basename, main

    # ------------------------------------------------------------------
    # exports
    # ------------------------------------------------------------------
    def export(self, output_path, fmt=None, **kwargs):
        '''
        Export with the format given by the extension of `output_path` (or by `fmt`).

        :param output_path: output path. ``.nc`` -> NetCDF file, ``.zarr`` -> Zarr store (with pyramids),
            no extension -> directory whose name is the base name of the files (NetCDF unless `fmt` says
            otherwise)
        :param fmt: ``'netcdf'`` or ``'zarr'`` to force the format (None: from the extension, else NetCDF).
            Must agree with the extension, if there is one.
        :param kwargs: passed to `export_to_netcdf` or `export_to_zarr`
        :return: for Zarr the path of the store, None for NetCDF
        '''
        from_name = _format_from_name(output_path)
        if fmt is not None and fmt not in _FORMAT_EXTENSION:
            raise ValueError(f"unknown format {fmt!r}, use 'netcdf' or 'zarr'")
        if fmt is not None and from_name is not None and fmt != from_name:
            raise ValueError(f"format {fmt!r} contradicts the extension of {output_path!r} ({from_name})")
        fmt = fmt or from_name or 'netcdf'

        if fmt == 'netcdf':
            return self.export_to_netcdf(output_path, **kwargs)
        return self.export_to_zarr(output_path, **kwargs)

    def export_to_netcdf(self, output_path, snap_compliant=False):
        '''
        Write the product as compressed, packed NetCDF files: the main product and the ancillary data
        (coarse resolution). Existing files are replaced.

        :param output_path: ``<dir>/<name>.nc``: main file ``<dir>/<name>.nc``, ancillary
            ``<dir>/<name>_anc.nc``. Without extension: directory whose name is the base name
            (``<output_path>/<basename>.nc`` and ``<output_path>/<basename>_anc.nc``, as before)
        :param snap_compliant: write ``Rrs`` as one variable per band so the file opens in SNAP
        '''
        logging.info('export into encoded netcdf')
        if _is_remote(output_path):
            raise ValueError('NetCDF can only be written to a local path; use a .zarr output for object storage')

        complevel = self.complevel
        ds, encoding = self._prepare_export(snap_compliant)
        for enc in encoding.values():
            enc.update({'zlib': True, 'complevel': complevel})

        # get file naming and create folder
        directory, basename, main_file = self._output_names(output_path, '.nc')
        anc_file = os.path.join(directory, basename + '_anc.nc')

        # clean up to avoid permission denied
        for path in (main_file, anc_file):
            if os.path.exists(path):
                os.remove(path)

        # export full raster data
        ds.to_netcdf(main_file, encoding=encoding, engine=NETCDF_ENGINE)

        # export ancillary data (coarse resolution)
        anc = self.ancillary.copy()
        self._scrub_attrs(anc)
        encoding = {}
        for variable in list(anc.keys()):
            encoding[variable] = {"zlib": True, "complevel": complevel, 'grid_mapping': 'spatial_ref'}

        anc.to_netcdf(anc_file, 'w',
                      encoding=encoding,
                      engine=NETCDF_ENGINE)

        return

    def export_to_zarr(self, output_path, snap_compliant=False, chunks=None, shards=None, sharding=True,
                       compression='zstd', clevel=5, storage_options=None, max_concurrency=None,
                       consolidate=True, pyramids=True, pyramid_levels=None,
                       pyramid_min_size=PYRAMID_MIN_SIZE):
        '''
        Write the product as one cloud-ready Zarr store with pyramids: the main product at full resolution
        in the group ``0``, the overviews (2x, 4x, ... coarser) in the groups ``1``, ``2``, ... and the
        ancillary data in the group ``ancillary``. The root group holds the product attributes and the
        ``multiscales`` description only. An existing store is replaced.

        Same variables, attributes and packing as `export_to_netcdf` (int16 + scale/offset, see
        `_base_encoding`). Cloud layout: sharded chunks, Blosc/zstd compression, consolidated metadata, direct
        write to object storage (see the module documentation). Overviews: mean for continuous variables,
        nearest for flags and masks.

        All the writes are evaluated together by a single ``dask.compute``, so the lazy processing graph
        runs once (the overviews are computed from the full resolution chunks while they are in memory) and
        its chunks are written in parallel.

        Parameters
        ----------
        output_path : str
            ``<dir>/<name>.zarr`` or a remote URL ``s3://bucket/<name>.zarr`` (the ``.zarr`` extension is
            required for a URL). Without extension the local store is ``<output_path>/<basename>.zarr``
            (as before).
        snap_compliant : bool, optional
            Write ``Rrs`` as one variable per band (kept for parity with NetCDF; SNAP itself reads NetCDF).
        chunks : dict, optional
            Stored chunk size per dimension, i.e. the unit of a read. Default 256 x 256 pixels with all
            wavelengths (512 x 512 when sharding is not available). Use ``'wl': 1`` to read single
            bands cheaply. Must divide the shard size, otherwise the whole shard is used.
        shards : dict, optional
            Size per dimension of the write blocks, which are the shards when sharding is on. Default
            ``{'y': 1024, 'x': 1024, 'wl': -1}``. Also the dask chunking of the computation: a shard is held
            in memory (float32, before packing) while it is computed, about ``1024 * 1024 * n_bands * 4``
            bytes per variable. Lower it (e.g. 512) when memory is tight.
        sharding : bool, optional
            Pack the chunks into shards (default True). Needs zarr-python >= 3 (Zarr format 3) and an xarray
            that accepts the ``shards`` encoding; otherwise a warning says what is missing and plain chunks
            are written (Zarr format 3 with zarr-python 3, format 2 with zarr-python 2).
        compression : str or None, optional
            Blosc codec: ``'zstd'`` (default), ``'lz4'`` (faster, larger), ``'zlib'``...; None keeps the
            Zarr default.
        clevel : int, optional
            Compression level (default 5).
        storage_options : dict, optional
            Credentials / options of the fsspec filesystem when `output_path` is a URL
            (e.g. ``{'profile': 'myprofile'}`` for S3, ``{'anon': False}``).
        max_concurrency : int, optional
            Number of concurrent I/O operations of zarr-python 3 (``async.concurrency``, default 10).
            Raising it (e.g. 32 or 64) speeds up the writes to object storage; ignored with zarr-python 2.
        consolidate : bool, optional
            Consolidate the metadata at the end (default True): the root group, which covers the whole
            hierarchy, and each level / ancillary group on its own. zarr-python 3 may warn that
            consolidated metadata is not part of the Zarr v3 specification; it is read transparently by
            ``xr.open_zarr`` and ignored by readers that do not support it.
        pyramids : bool, optional
            Write the overviews (default True).
        pyramid_levels : int, optional
            Number of overviews (default: as many as needed to get a coarsest level of at most
            `pyramid_min_size` pixels on its longest side).
        pyramid_min_size : int, optional
            See `pyramid_levels` (default 512).

        Returns
        -------
        str
            Path / URL of the Zarr store.
        '''
        import dask

        logging.info('export into encoded zarr')
        zarr_major = _zarr_major()
        zarr3 = (zarr_major or 0) >= 3  # Zarr format 3 (default of zarr-python 3)
        sharded = bool(sharding) and zarr3 and _xarray_supports_shards()
        if sharding and not sharded:
            reason = 'zarr-python >= 3' if not zarr3 else 'a recent xarray (the installed one has no `shards` encoding)'
            logging.warning(f'sharding needs {reason}: writing plain chunks')
        chunks = dict(chunks) if chunks else dict(_DEFAULT_CHUNKS_SHARDED if sharded else _DEFAULT_CHUNKS)
        block = dict(shards) if shards else dict(_DEFAULT_SHARDS)
        compression_enc = _compression_encoding(compression, clevel, zarr_major)

        ds, encoding = self._prepare_export(snap_compliant)
        # encodings inherited from the input (e.g. NetCDF or Zarr chunks, codecs) must not leak into the output,
        # including those of the coordinates (time, spatial_ref)
        _reset_encoding(ds)
        ds = _load_coords(ds.chunk({d: s for d, s in block.items() if d in ds.dims}))

        # lazy pyramid: each level is the coarsening of the previous one
        levels = [ds]
        if pyramids:
            n_over = _pyramid_levels(ds.sizes['y'], ds.sizes['x'], pyramid_min_size, pyramid_levels)
            for _ in range(n_over):
                level = _coarsen(levels[-1])
                level = _load_coords(level.chunk({d: s for d, s in block.items() if d in level.dims}))
                levels.append(level)

        # description of the levels, stored in the root attributes
        groups = [LEVEL_GROUP.format(k) for k in range(len(levels))]
        layout = [{'group': group,
                   'scale': PYRAMID_FACTOR ** k,
                   'shape': [int(lv.sizes['y']), int(lv.sizes['x'])]}
                  for k, (group, lv) in enumerate(zip(groups, levels))]
        # root group: no data, product attributes + description of the levels for viewers
        root = xr.Dataset(attrs=dict(ds.attrs))
        root.attrs['multiscales'] = [{
            'description': 'pyramid of the main product, group k is {:d}**k coarser than group 0'.format(PYRAMID_FACTOR),
            'resampling': 'mean ignoring NaN for continuous variables, nearest for integer variables (flags, mask)',
            'layout': layout}]

        def level_encoding(level):
            '''Encoding of one level: packing + compression + stored chunks (+ shards = dask blocks).'''
            enc = {var: dict(e) for var, e in _coord_encoding().items() if var in level.variables}
            for var, e in encoding.items():
                if var not in level.variables:
                    continue
                e = dict(e)
                e.update(compression_enc)
                if level[var].dims and level[var].chunks is not None:
                    inner, blk = _var_chunking(level[var].dims, level[var].chunks, chunks)
                    e['chunks'] = inner
                    if sharded and inner != blk:
                        e['shards'] = blk
                enc[var] = e
            return enc

        anc = self.ancillary.copy()
        self._scrub_attrs(anc)
        _reset_encoding(anc)
        anc = _load_coords(anc.chunk(-1))  # coarse data: one chunk per variable
        encoding_anc = {var: dict(e) for var, e in _coord_encoding().items() if var in anc.variables}
        encoding_anc.update({var: dict(compression_enc, grid_mapping='spatial_ref',
                                       chunks=tuple(c[0] for c in anc[var].chunks))
                             for var in anc.data_vars})

        # target: URL (written directly to the bucket) or local path
        if _is_remote(output_path):
            store = str(output_path).rstrip('/')
            if _format_from_name(store) != 'zarr':
                raise ValueError(f"the URL of a Zarr store has to end with '.zarr': {output_path!r}")
        else:
            directory, basename, store = self._output_names(output_path, '.zarr')
        _clear_store(store, storage_options)

        write_kwargs = {'compute': False, 'consolidated': False}
        if zarr3:
            write_kwargs['zarr_format'] = 3  # explicit: sharding is a Zarr format 3 feature
        if storage_options:
            write_kwargs['storage_options'] = storage_options

        # metadata is written now, the data when the delayed objects are computed
        delayed = [root.to_zarr(store, mode='w', **write_kwargs)]
        for group, level in zip(groups, levels):
            delayed.append(level.to_zarr(store, group=group, mode='a',
                                         encoding=level_encoding(level), **write_kwargs))
        delayed.append(anc.to_zarr(store, group=ANCILLARY_GROUP, mode='a', encoding=encoding_anc, **write_kwargs))
        # one evaluation of the whole graph; zarr's own I/O concurrency can be raised for object storage
        if max_concurrency and zarr3:
            import zarr
            with zarr.config.set({'async.concurrency': int(max_concurrency)}):
                dask.compute(*delayed)
        else:
            dask.compute(*delayed)

        if consolidate:
            try:
                _consolidate(store, storage_options, zarr_major, groups + [ANCILLARY_GROUP])
            except Exception as e:
                logging.warning(f'could not consolidate zarr metadata: {e}')

        return store