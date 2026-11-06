Usage
=====

GRS needs, for each image:

* one Level-1 product: Sentinel-2 ``MSIL1C`` (``.SAFE``), Landsat 8/9 Level-1, or a Level-1 product
  already converted to NetCDF by `GRSdriver <https://github.com/Tristanovsk/GRSdriver>`__;
* the CAMS atmospheric composition file covering the date and the area (see the
  :doc:`CAMS tutorial <examples/auxdata>`);
* optionally a DEM and a SURFWATER water mask, already subset to the tile.

Command line
------------

.. literalinclude:: ../../grs/run.py
   :language: text
   :start-at: Usage:
   :end-before: '''

The command line writes the product in ``<odir>/<name>/`` with ``<name>`` the input basename where
``L1C``/``L1GT``/``L1TP`` is replaced by ``--levname`` (default ``L2AGRS``), without extension, plus
``--suffix`` (default ``_V<GRS version>``). Sentinel-2 images with a cloud cover above
``--max_cloud_cover`` are skipped.

Python
------

.. code-block:: python

   import grs

   process = grs.Process()
   process.execute('S2A_MSIL1C_20181019T102031_N0500_R065_T30PYT_20230815T043850.SAFE',
                   cams_file='cams_forecast_2018-10.nc',
                   dem_file='COP-DEM_GLO-30-DGED_30PYT.tif',
                   resolution=20,
                   odir='S2A_MSIL2AGRS_20181019T102031_N0500_R065_T30PYT_20230815T043850')

   # process.l2a.l2_prod is the Level-2 xarray.Dataset, held in memory
   process.write_output()

The options of :py:meth:`~grs.grs_process.Process.execute`:

``resolution``
   Pixel size of the processing (10, 20 or 60 m for Sentinel-2).
``cams_file``
   CAMS file of the date (surface pressure, gas columns, aerosol optical depth, wind).
``dem_file``
   DEM (GeoTIFF, subset to the tile) used to compute the surface pressure; without it, the surface
   pressure of CAMS is used.
``surfwater_file``
   SURFWATER water mask (GeoTIFF) used for the pixel classification.
``opac_model``
   Force the aerosol model (``'ARCT_rh70'``, ``'COAV_rh70'``, ``'DESE_rh70'``, ``'MACL_rh70'``,
   ``'URBA_rh70'``); by default the model closest to the spectral AOD of CAMS is used.
``scale_aot``
   Scaling factor of the CAMS AOD at 550 nm.
``allpixels``
   Process all the pixels, not only the water pixels.
``snap_compliant``
   Write one ``Rrs_<wl>`` variable per band, readable by ESA SNAP.

Output product
--------------

:py:meth:`~grs.grs_process.Process.write_output` writes two NetCDF files in the output directory
``odir``: ``<name>.nc`` (main product, full resolution) and ``<name>_anc.nc`` (ancillary data on the
CAMS grid, dimensions ``xc``, ``yc``), with ``<name>`` the name of the directory. See
:doc:`api/output` for the layout.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Variable
     - Dimensions
     - Description
   * - ``Rrs``
     - wl, y, x
     - remote-sensing reflectance (sr\ :sup:`-1`)
   * - ``BRDFg``
     - y, x
     - sunglint reflectance in the SWIR (see :doc:`methods`)
   * - ``aot550``
     - y, x
     - aerosol optical thickness at 550 nm (from CAMS)
   * - ``flags``
     - y, x
     - pixel classification bit field (see :doc:`api/mask`)
   * - ``mask``
     - y, x
     - 0 for the good water pixels, 1 elsewhere
   * - ``sza``, ``vza``, ``raa``
     - y, x
     - solar zenith, viewing zenith and relative azimuth angles (degrees)
   * - ``dem``
     - y, x
     - elevation (m), empty without ``dem_file``

The global attributes record the metadata of the Level-1 product, the masking thresholds and the
statistics of the flags.

Read the product back with:

.. code-block:: python

   import xarray as xr

   name = 'S2A_MSIL2AGRS_20181019T102031_N0500_R065_T30PYT_20230815T043850'
   l2a = xr.open_dataset(f'{name}/{name}.nc')
   anc = xr.open_dataset(f'{name}/{name}_anc.nc')
