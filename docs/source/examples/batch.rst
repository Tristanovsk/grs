Batch processing
================

The notebook below processes a time series of Sentinel-2 L1C images (or Landsat Collection 2 Level-1
images) of one or several tiles into GRS L2A products, with
:py:class:`grs.Process <grs.grs_process.Process>`:

- configuration of the batch: input and output folders, tiles, date range, auxiliary data and
  processing options, and check of the look-up tables set in ``grs/config.yml``;
- list of the Level-1 images found on disk (``<idir>/<tile>/YYYY/MM/DD/<image>``) with their cloud
  cover, and selection of a date range;
- CAMS file of the month of each image (see the :doc:`CAMS tutorial <auxdata>`), optional DEM of the
  tile, and output names as with the ``grs`` command line;
- the batch itself: images already processed, too cloudy or without CAMS file are skipped, a failing
  image is logged and does not stop the batch, and a summary table is saved;
- a quick look at one L2A product (NetCDF);
- the same batch from the command line, one ``grs`` process per image, which is more robust for long
  batches.

.. code-block:: python

   import grs

   for date, row in selection.iterrows():          # one row per Level-1 image
       process = grs.Process()
       process.execute(row['path'], odir=row['outdir'], cams_file=row['cams_file'],
                       dem_file=row['dem_file'], resolution=20)
       process.write_output()                      # <outdir>/<name>.nc and <name>_anc.nc

The notebook is rendered with the outputs of its last run: the images of 1 to 4 August 2025 of tile
31TEJ (Gulf of Lion: Thau lagoon and the coast of Montpellier). With GRS 2.2.0, a full tile at 20 m
takes about 8 minutes and up to 30 GB of memory. The input images are not distributed with GRS: set the
paths of the configuration cell to run it on your own data.

.. toctree::
   :maxdepth: 1
   :caption: Batch processing

   batch/grs_s2_batch
