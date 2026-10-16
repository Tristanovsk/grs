Usage
=====

GRS needs, for each image:

* one Level-1 product: Sentinel-2 ``MSIL1C`` (``.SAFE``) or Landsat 8/9 Level-1;
* the CAMS atmospheric composition file covering the date and the area (see the
  :doc:`CAMS tutorial <examples/auxdata>`);
* optionally a DEM and a SURFWATER water mask, already subset to the tile.

Command line
------------

.. literalinclude:: ../../grs/run.py
   :language: text
   :start-at: Usage:
   :end-before: '''

To process many images, see the :doc:`batch processing tutorial <examples/batch>`: it lists the
Level-1 images of a date range and runs them from Python, or writes one ``grs`` command per image to
run them from a shell (several at a time if the memory allows it).

Python
------

.. code-block:: python

   import grs

   process = grs.Process()
   process.execute('S2A_MSIL1C_20181019T102031_N0500_R065_T30PYT_20230815T043850.SAFE',
                   cams_file='cams_forecast_2018-10.nc',
                   dem_file='COP-DEM_GLO-30-DGED_30PYT.tif',
                   resolution=20,
                   aot_resolution=240,
                   odir='S2A_MSIL1C_20181019T102031_L2AGRS.zarr')

   # the processing is lazy: process.l2a.l2_prod is a dask-backed xarray.Dataset,
   # computed chunk by chunk when it is written
   process.write_output()

The main options of :py:meth:`~grs.grs_process.Process.execute`:

``resolution``
   Pixel size of the processing (10, 20 or 60 m for Sentinel-2).
``aot_resolution``
   Pixel size of the coarse grid of the aerosol retrieval (default 240 m); must be a multiple of
   ``resolution`` (e.g. 60, 120, 240, 480 m). See :ref:`aerosol retrieval <methods-aerosol>`.
``aot_chunk``, ``aot_smoothing_blocks``, ``min_pixel_prop``, ``outlier_factor``, ``fine_step``, ``fine_window``, ``aot_step``
   The other settings of the aerosol retrieval. Their defaults, like the one of ``aot_resolution``, are set
   in the ``aerosol`` section of the configuration (see :doc:`configuration`); the arguments of
   :py:meth:`~grs.grs_process.Process.execute` override them.
``opac_model``
   Force the aerosol model (``'ARCT_rh70'``, ``'COAV_rh70'``, ``'DESE_rh70'``, ``'MACL_rh70'``,
   ``'URBA_rh70'``); by default the model closest to the spectral AOD of CAMS is used.
``scale_aot``
   Scaling factor of the CAMS AOD.
``allpixels``
   Process all the pixels, not only the water pixels.
``chunk``
   Block size (pixels) of the final processing; a 512 x 512 block needs about 350 MB per thread.

Output product
--------------

The format is given by the extension of the output path: ``.zarr`` (default, with overviews, can be
written to object storage) or ``.nc`` (NetCDF, with a separate ``_anc.nc`` file for the ancillary data).
See :doc:`api/output` for the layouts.

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
     - aerosol optical thickness at 550 nm
   * - ``wind``
     - y, x
     - wind speed used for the sunglint (m s\ :sup:`-1`)
   * - ``flags``
     - y, x
     - pixel classification bit field (see :doc:`api/mask`)
   * - ``mask``
     - y, x
     - 0 for the good water pixels, 1 elsewhere
   * - ``sza``, ``vza``, ``raa``
     - y, x
     - solar zenith, viewing zenith and relative azimuth angles (degrees)

The global attributes record the processing settings, among them ``opac_model`` and
``aot_resolution``.

Read a Zarr product back with:

.. code-block:: python

   import xarray as xr

   l2a = xr.open_zarr('S2A_MSIL1C_20181019T102031_L2AGRS.zarr', group='0')        # full resolution
   anc = xr.open_zarr('S2A_MSIL1C_20181019T102031_L2AGRS.zarr', group='ancillary')
