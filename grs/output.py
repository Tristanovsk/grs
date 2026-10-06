'''
Module to set output product and format (L2A): NetCDF and Zarr export.

An :class:`L2aProduct` gathers the L2A variables (remote-sensing reflectance, aerosol, geometry, flags,
DEM, ...) and the coarse-resolution ancillary data (CAMS fields, gaseous transmittance), and writes them
with a fixed, lossy-but-controlled packing (``int16`` with scale/offset, see `_base_encoding`).

Output layouts
==============

``export_to_netcdf(output_path)``  (unchanged layout)::

    <output_path>/<basename>.nc        main product
    <output_path>/<basename>_anc.nc    ancillary data (coarse grid, dimensions xc, yc)

``export_to_zarr(output_path)``::

    <output_path>/<basename>.zarr/            root group = main product
    <output_path>/<basename>.zarr/ancillary/  group = ancillary data (coarse grid, dimensions xc, yc)

Read back with ``xr.open_zarr(store)`` (main) and ``xr.open_zarr(store, group='ancillary')``.
The packed variables are unpacked automatically (``scale_factor`` / ``add_offset`` / ``_FillValue``).

Zarr export and dask
====================

The product is normally dask-backed (lazy). The Zarr export is built to evaluate it **once**:

* the main and ancillary datasets are rechunked to regular chunks (default ``y = x = 1024`` and all
  wavelengths in one chunk, the same layout ``Product`` reads) so Zarr chunks and dask chunks coincide;
* both writes are prepared with ``compute=False`` and run together in one ``dask.compute`` call, so
  intermediate results shared by the main and the ancillary graphs (e.g. the gaseous transmittance)
  are computed only once and chunks are written in parallel as soon as they are ready;
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
        (as strings), the geometry angles (``vza``, ``raa`` averaged over wavelengths if band dependent;
        ``sza``), ``flags`` and ``dem``; spatial dimensions ``x``/``y``, CRS written.
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

        self.ancillary = xr.merge([transmittance_raster, cams_raster], compat='override')
        self.ancillary = self.ancillary.rio.set_spatial_dims(x_dim='xc', y_dim='yc')
        self.ancillary = self.ancillary.rio.write_coordinate_system()
        self.ancillary = self.ancillary.rio.write_crs()

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
    def _output_names(output_path):
        '''Create the output directory; return ``(directory, basename)`` (robust to a trailing slash).'''
        output_path = os.path.normpath(output_path)
        os.makedirs(output_path, exist_ok=True)
        return output_path, os.path.basename(output_path)

    # ------------------------------------------------------------------
    # exports
    # ------------------------------------------------------------------
    def export(self, output_path, fmt='netcdf', **kwargs):
        '''
        Export with the chosen format.

        :param output_path: output directory (its name is used as base name of the files)
        :param fmt: ``'netcdf'`` or ``'zarr'``
        :param kwargs: passed to `export_to_netcdf` or `export_to_zarr`
        '''
        if fmt == 'netcdf':
            return self.export_to_netcdf(output_path, **kwargs)
        if fmt == 'zarr':
            return self.export_to_zarr(output_path, **kwargs)
        raise ValueError(f"unknown format {fmt!r}, use 'netcdf' or 'zarr'")

    def export_to_netcdf(self, output_path, snap_compliant=False):
        '''
        Write the product as compressed, packed NetCDF files: ``<output_path>/<basename>.nc`` (main) and
        ``<output_path>/<basename>_anc.nc`` (ancillary, coarse resolution). Existing files are replaced.

        :param output_path: output directory (its name is used as base name of the files)
        :param snap_compliant: write ``Rrs`` as one variable per band so the file opens in SNAP
        '''
        logging.info('export into encoded netcdf')

        complevel = self.complevel
        ds, encoding = self._prepare_export(snap_compliant)
        for enc in encoding.values():
            enc.update({'zlib': True, 'complevel': complevel})

        # get file naming and create folder
        output_path, basename = self._output_names(output_path)
        ofile = os.path.join(output_path, basename)

        # clean up to avoid permission denied
        if os.path.exists(ofile + '.nc'):
            os.remove(ofile + '.nc')
        if os.path.exists(ofile + '_anc.nc'):
            os.remove(ofile + '_anc.nc')

        # export full raster data
        ds.to_netcdf(ofile + '.nc', encoding=encoding, engine=NETCDF_ENGINE)

        # export ancillary data (coarse resolution)
        anc = self.ancillary.copy()
        self._scrub_attrs(anc)
        encoding = {}
        for variable in list(anc.keys()):
            encoding[variable] = {"zlib": True, "complevel": complevel, 'grid_mapping': 'spatial_ref'}

        anc.to_netcdf(ofile + '_anc.nc', 'w',
                      encoding=encoding,
                      engine=NETCDF_ENGINE)

        return

    def export_to_zarr(self, output_path, snap_compliant=False, chunks=None, consolidate=True):
        '''
        Write the product as one Zarr store: ``<output_path>/<basename>.zarr`` with the main product in the
        root group and the ancillary data in the group ``ancillary``. An existing store is replaced.

        Same variables, attributes and packing as `export_to_netcdf` (int16 + scale/offset, see
        `_base_encoding`); the compression is Zarr's default compressor (the codec option names differ
        between zarr-python 2 and 3, so none is forced).

        The main and the ancillary writes are evaluated together by a single ``dask.compute``, so the lazy
        processing graph runs once and its chunks are written in parallel.

        Parameters
        ----------
        output_path : str
            Output directory (its name is used as base name of the store).
        snap_compliant : bool, optional
            Write ``Rrs`` as one variable per band (kept for parity with NetCDF; SNAP itself reads NetCDF).
        chunks : dict, optional
            Chunk sizes of the main product per dimension. Default ``{'y': 1024, 'x': 1024, 'wl': -1}``
            (all bands in one chunk). Use smaller spatial chunks (e.g. 512) for web viewing, ``'wl': 1`` to
            read single bands cheaply. The ancillary data is written in one chunk per variable.
        consolidate : bool, optional
            Consolidate the metadata of the whole store at the end (default True).

        Returns
        -------
        str
            Path of the Zarr store.
        '''
        import dask

        logging.info('export into encoded zarr')
        if chunks is None:
            chunks = {'y': 1024, 'x': 1024, 'wl': -1}

        ds, encoding = self._prepare_export(snap_compliant)
        # variable encodings inherited from the input (e.g. NetCDF or Zarr chunks) must not leak into the output
        for var in ds.data_vars:
            ds[var].encoding = {}
        ds = ds.chunk({d: s for d, s in chunks.items() if d in ds.dims})
        for var, enc in encoding.items():
            if ds[var].chunks is not None:
                enc['chunks'] = tuple(c[0] for c in ds[var].chunks)

        anc = self.ancillary.copy()
        self._scrub_attrs(anc)
        for var in anc.data_vars:
            anc[var].encoding = {}
        anc = anc.chunk(-1)
        encoding_anc = {var: {'grid_mapping': 'spatial_ref',
                              'chunks': tuple(c[0] for c in anc[var].chunks)}
                        for var in anc.data_vars}

        output_path, basename = self._output_names(output_path)
        store = os.path.join(output_path, basename + '.zarr')
        if os.path.exists(store):
            shutil.rmtree(store)

        # metadata is written now, the data when the delayed objects are computed
        delayed_main = ds.to_zarr(store, mode='w', encoding=encoding,
                                  compute=False, consolidated=False)
        delayed_anc = anc.to_zarr(store, group='ancillary', mode='a', encoding=encoding_anc,
                                  compute=False, consolidated=False)
        dask.compute(delayed_main, delayed_anc)  # one evaluation of the whole graph

        if consolidate:
            try:
                import zarr
                zarr.consolidate_metadata(store)
            except Exception as e:
                logging.warning(f'could not consolidate zarr metadata: {e}')

        return store