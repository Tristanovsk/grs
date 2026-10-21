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

The output granule is a single Level-2 GRS product covering the same footprint and resolution as
the input tile: a Zarr store with pyramids (default), or NetCDF files (see `Output product format`_).

Scheduling and Triggers
~~~~~~~~~~~~~~~~~~~~~~~~

GRSProcessor has no built-in scheduler; it is triggered externally, granule by granule:

- **Interactively / single granule**: manual call to ``grs <input_file> ...`` (see
  :doc:`usage` for CLI examples).
- **Batch, on an HPC cluster (SLURM)**: submitted as a job via ``sbatch``, see
  `ci-func-run.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-func-run.slurm>`_
  (runs ``grs`` / ``exe/launcher.py`` directly on one or more tiles) and
  `ci-docker-run.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-docker-run.slurm>`_
  (same, inside the Singularity/Docker container).
  The CI/validation pipeline (see ``.gitlab-ci.yml``) chains jobs with
  ``sbatch --depend=afterany:<jobid>``:
  `ci-init_env.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-init_env.slurm>`_
  (environment setup) →
  `ci-func-run.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-func-run.slurm>`_ →
  `ci-cleanup_env.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-cleanup_env.slurm>`_.
  Auxiliary CAMS retrieval also runs on SLURM
  (`download_cams.slurm <https://github.com/CNES/GRSprocessor/blob/main/ecmwf/download_cams.slurm>`_,
  `slurm_cds.slurm <https://github.com/CNES/GRSprocessor/blob/main/ecmwf/slurm_cds.slurm>`_).
  There is no periodicity of its own; scheduling (e.g. reprocessing on new L1C acquisitions)
  is delegated to the job scheduler / calling scripts.

  .. note::
     Legacy PBS launchers (``grs_launcher.pbs``, ``grs_mpi_launcher.pbs``,
     ``launch_grs_exemple.pbs``) are still present at the repository root but are
     superseded by SLURM on the current CNES cluster.
- **Containerized**: via the Docker image built from
  `Dockerfile <https://github.com/CNES/GRSprocessor/blob/main/Dockerfile>`_
  (see :doc:`index` "Compile Docker image locally").

Dataflow
--------

.. only:: not latex

   .. mermaid:: _diagrams/dataflow.mmd

.. only:: latex

   The dataflow diagram is drawn with Mermaid: see the
   `online documentation <https://grs.readthedocs.io/en/v3.x/processing_chain.html#dataflow>`__.

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
   * - Level-2 GRS product (Zarr)
     - 1..1
     - Water-leaving reflectance product, one Zarr store per granule with the full resolution, its
       overviews (pyramids) and the ancillary data (see `Output product format`_ and
       `Data Types`_ for naming). Default output of the ``grs`` command line.
   * - Level-2 GRS product (NetCDF)
     - 0..1
     - Alternative format, from the Python API (output path ending in ``.nc``): main file and
       ancillary file ``<name>_anc.nc``, without overviews.
   * - Run log (``log_file.log``)
     - 1..1
     - Full execution log (all levels), written next to the output product.
   * - Error log (``error.log``)
     - 0..1
     - Subset of the run log containing ``ERROR``-level records only; empty if the run
       completed without error.

Output product format
~~~~~~~~~~~~~~~~~~~~~

The format is chosen from the output path given to
:py:meth:`Process.write_output <grs.grs_process.Process.write_output>`: ``.zarr`` gives a Zarr store,
``.nc`` NetCDF files; a path without extension (as used by the ``grs`` command line) is a directory
that receives a Zarr store named after it.

**Zarr store with pyramids** (default)::

    <name>.zarr/              root group: product attributes and ``multiscales``, no data
    <name>.zarr/0/            full resolution (level 0)
    <name>.zarr/1/            overview, 2 x coarser
    <name>.zarr/2/            overview, 4 x coarser
    ...
    <name>.zarr/ancillary/    ancillary data on the coarse grid (dimensions xc, yc)

* **Pyramids**: each level is built from the previous one by 2 x 2 blocks: mean (ignoring NaN) for the
  continuous variables (``Rrs``, ``aot550``, angles...), top-left pixel for ``flags`` and ``mask``
  (averaging bit flags would give meaningless values: use the full resolution for any statistics on
  the flags). Levels are added until the longest side is at most 512 pixels; they add about 1/3 of the
  full-resolution size. The root attribute ``multiscales`` lists the levels (group, scale, shape) for
  viewers, and the georeferencing is updated for each level.
* **Cloud-ready layout**: inner chunks of 256 x 256 pixels (all wavelengths) packed into 1024 x 1024
  shards (Zarr format 3, zarr-python >= 3; plain 512 x 512 chunks otherwise), Blosc/zstd compression and
  consolidated metadata, so that a window of the image can be read without downloading the rest.
* **Remote output**: the store can be written directly to object storage (``s3://``, ``gs://``,
  ``az://`` URLs with the matching fsspec driver).
* **Packing**: the variables are stored as scaled ``int16`` (e.g. ``Rrs``: scale 1e-5, offset 0.3),
  unpacked automatically when read.
* **Single pass**: the processing is lazy (dask) and the full resolution, the overviews and the
  ancillary data are written in one computation, chunk by chunk, so the full image is never in memory.

Read the product back with xarray:

.. code-block:: python

   import xarray as xr

   l2a = xr.open_zarr('<name>.zarr', group='0')          # full resolution
   quicklook = xr.open_zarr('<name>.zarr', group='3')    # 8 x coarser
   anc = xr.open_zarr('<name>.zarr', group='ancillary')

**NetCDF**: ``<name>.nc`` (main product) and ``<name>_anc.nc`` (ancillary data), local path only, no
overviews; ``snap_compliant=True`` writes one variable per band for ESA SNAP.

See :doc:`api/output` for all the export options (chunk and shard sizes, compression, number of
levels...).

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
the SLURM job scripts shipped in the repository:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Resource Type
     - Quantity
     - Source / notes
   * - CPU
     - 8 tasks (``-n 8``, single node) for a direct/native run — 4 tasks for the
       containerized (Singularity) run
     - `ci-func-run.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-func-run.slurm>`_,
       `ci-docker-run.slurm <https://github.com/CNES/GRSprocessor/blob/main/grs/tests/ci/ci-docker-run.slurm>`_
   * - RAM
     - 4 GB/task (32 GB total) for the direct run — 32 GB/task (128 GB total) for the
       containerized run
     - same SLURM scripts (``--mem-per-cpu``)
   * - Execution time
     - walltime budget: 30 min (direct run) — 2 h (containerized run)
     - same SLURM scripts (``--time``); actual runtime per granule is much shorter and
       logged as ``total_run_time`` (see `Log Format`_)
   * - Disk — input
     - size of one L1C/L1 product (SAFE or tar) plus the CAMS netCDF file
     - depends on sensor/resolution
   * - Disk — output
     - one Level-2 Zarr store (full resolution + about 1/3 for the overviews) + ``log_file.log`` +
       ``error.log`` per granule
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
     - Zarr store (with pyramids) of the water-leaving reflectance product
     - one tile, one date
     - ``<odir>/<name>/<name>.zarr`` with ``<name>`` = input basename with L1C/L1GT/L1TP replaced
       by ``levname`` (default L2AGRS), extension removed, plus ``suffix`` (default
       ``_V<GRS version>``), see ``exe.procutils.misc.set_ofile``. The logs are written in
       ``<odir>/<name>/``.
