Processing Chain Description
=============================

This page documents GRSProcessor following the
`Processing Chain Documentation Template <https://processing-chain-guidelines.readthedocs.io/en/latest/Data_processing_chain_template/>`_.

Process description
--------------------

Description
~~~~~~~~~~~

GRS (Glint Removal for Sentinel-2-like sensors) is an atmospheric and sunglint correction
processor for high-spatial-resolution, multispectral optical satellite images acquired over
aquatic environments. Given one Level-1 (top-of-atmosphere) product, it corrects successively
for gaseous absorption, diffuse sky light reflected by the air-water interface, and sunglint,
to retrieve the water-leaving radiance/reflectance at the surface level. The full scientific
description is available on the :doc:`index` page.

The processor is invoked through the ``grs`` command-line executable
(entry point ``grs.run:main``, see :py:mod:`grs.run`), which wraps the core
:py:class:`grs.grs_process.Process` class.

Application Domain (Granule)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One processing run (one granule) corresponds to **one Level-1 product for one acquisition
date and one tile/footprint**:

- one Sentinel-2 A/B ``MSIL1C`` ``.SAFE`` product, or
- one Landsat-8/9 Level-1 (``L1TP``/``L1GT``) product.

The output granule is a single Level-2 GRS product (netCDF) covering the same footprint and
resolution as the input tile.

Scheduling and Triggers
~~~~~~~~~~~~~~~~~~~~~~~~

GRSProcessor has no built-in scheduler; it is triggered externally, granule by granule:

- **Interactively / single granule**: manual call to ``grs <input_file> ...`` (see
  :doc:`index` "Testing" section for CLI examples).
- **Batch, on an HPC cluster (PBS)**: submitted as a job via ``qsub``, see
  `grs_launcher.pbs <https://github.com/CNES/GRSprocessor/blob/main/grs_launcher.pbs>`_
  (single tile) and
  `grs_mpi_launcher.pbs <https://github.com/CNES/GRSprocessor/blob/main/grs_mpi_launcher.pbs>`_
  (list of tiles/dates processed in parallel, driven by a CSV list file).
  There is no periodicity of its own; scheduling (e.g. reprocessing on new L1C acquisitions)
  is delegated to the job scheduler / calling scripts.
- **Containerized**: via the Docker image built from
  `Dockerfile <https://github.com/CNES/GRSprocessor/blob/main/Dockerfile>`_
  (see :doc:`index` "Compile Docker image locally").

Inputs
------

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Data Type Name
     - Cardinality
     - Selection Criteria
   * - Level-1 product (``<input_file>``)
     - 1..1 (mandatory)
     - One Sentinel-2 A/B ``MSIL1C`` ``.SAFE`` directory, or one Landsat-8/9 L1 product
       (e.g. a ``.tar``), passed as the positional CLI argument.
   * - CAMS file (``--cams_file``)
     - 1..1 (mandatory)
     - Copernicus Atmosphere Monitoring Service netCDF file covering the acquisition date/area,
       used to interpolate atmospheric pressure and aerosol optical thickness (AOT).
   * - DEM file (``--dem_file``)
     - 0..1 (optional)
     - GeoTIFF digital elevation model already subset to the input tile footprint; used to
       rescale the Rayleigh optical thickness with altitude.
   * - Surface water mask (``--surfwater``)
     - 0..1 (optional)
     - GeoTIFF surface-water mask (e.g. SURFWATER product) restricting/flagging the processed
       water pixels.

Outputs
-------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Data Type Name
     - Cardinality
     - Description
   * - Level-2 GRS product (netCDF)
     - 1..1
     - Water-leaving radiance/reflectance product, one file per granule (see
       `Data Types`_ below for naming).
   * - Run log (``log_file.log``)
     - 1..1
     - Full execution log (all levels), written next to the output product.
   * - Error log (``error.log``)
     - 0..1
     - Subset of the run log containing ``ERROR``-level records only; empty if the run
       completed without error.

Return Codes
------------

GRSProcessor is a Python CLI (``docopt``-based); process errors are primarily reported through
the logs rather than through a rich set of process exit codes. Currently observed behavior of
``grs`` (see :py:func:`grs.run.main`):

- ``0``: normal process exit. This includes the case where processing raised an exception —
  the exception is caught, logged as ``ERROR`` to ``error.log`` (with traceback), and the
  process exits normally. **Callers must inspect ``error.log`` for emptiness rather than rely
  on the exit code to detect failure.**
- ``-1``: the granule was skipped because the output product already exists and
  ``--no_clobber`` was set.
- Non-zero interpreter exit codes may also occur for errors raised before logging is
  initialized (e.g. invalid CLI arguments rejected by ``docopt``, unreadable input product).

Log Format
----------

Logging is configured by :py:class:`grs.class_logger.ServiceLogger`, one log file and one
error-only log file per run (see `Outputs`_).

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Description
   * - ``asctime``
     - Timestamp, ``YYYY-MM-DDThh:mm:ss.mmm``.
   * - ``levelname``
     - Classification level: ``DEBUG``, ``INFO``, ``WARNING`` or ``ERROR``.
   * - ``name``
     - Logger/class name (e.g. module emitting the record).
   * - ``funcName``
     - Method/function name that emitted the record.
   * - ``message``
     - Free-text message body.

Format string: ``"%(asctime)s.%(msecs)03d  %(levelname)s | %(name)s::%(funcName)s | %(message)s"``.

At the end of a run, GRSProcessor also logs resource-usage statistics (see
`Required Resources`_) as additional ``INFO`` records: ``max_rss``, ``sys_cpu``, ``user_cpu``,
``total_run_time``, and per-level message counters (``error``, ``warning``, ``info``, ``debug``).

Required Resources
-------------------

Actual needs depend on product resolution (10/20/60 m) and tile size. Indicative values from
the HPC job templates shipped in the repository:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Resource Type
     - Quantity
     - Source / notes
   * - CPU
     - 16 cores (single tile) — up to 40 cores / 8 MPI processes (batch list of tiles)
     - `grs_launcher.pbs <https://github.com/CNES/GRSprocessor/blob/main/grs_launcher.pbs>`_,
       `grs_mpi_launcher.pbs <https://github.com/CNES/GRSprocessor/blob/main/grs_mpi_launcher.pbs>`_
   * - RAM
     - 32 GB (single tile) — up to 180 GB (batch)
     - same PBS templates
   * - Execution time
     - walltime budget: 24 h (single tile) — 48 h (batch)
     - same PBS templates; actual runtime per granule is much shorter and logged as
       ``total_run_time`` (see `Log Format`_)
   * - Disk — input
     - size of one L1C/L1 product (SAFE or tar) plus the CAMS netCDF file
     - depends on sensor/resolution
   * - Disk — output
     - one netCDF Level-2 product + ``log_file.log`` + ``error.log`` per granule
     - written under ``--odir``

Data Types
----------

.. list-table::
   :header-rows: 1
   :widths: 18 32 15 35

   * - Data Name
     - Description
     - Granule
     - Nomenclature
   * - Level-1 input product
     - Sentinel-2 or Landsat-8/9 top-of-atmosphere product
     - one tile, one date
     - Sentinel-2: ``S2{A,B}_MSIL1C_<datetime>_N<baseline>_R<orbit>_T<tile>_<datetime>.SAFE``.
       Landsat: standard USGS Collection-2 L1 product naming (``L1TP``/``L1GT``).
   * - CAMS file
     - Copernicus Atmosphere Monitoring Service atmospheric composition forecast
     - regional/global grid, one date
     - e.g. ``YYYY-MM-DD-cams-global-atmospheric-composition-forecasts.nc``
   * - DEM file
     - GeoTIFF elevation raster subset to the tile footprint
     - one tile
     - project/site-specific, e.g. ``COP-DEM_GLO-30-DGED_<tile>.tif``
   * - Surface water mask
     - GeoTIFF water/non-water mask
     - one tile, one date
     - e.g. SURFWATER product naming,
       ``SURFWATER_OPTICAL-SINGLE_<tile>_<datetime>_<datetime>_<version>.tif``
   * - Level-2 GRS output
     - netCDF water-leaving radiance/reflectance product
     - one tile, one date
     - ``<basename with L1C/L1GT/L1TP replaced by levname (default L2AGRS)><suffix>.nc``,
       default ``suffix`` is ``_V<GRS version>`` (see
       :py:func:`exe.procutils.misc.set_ofile`)
