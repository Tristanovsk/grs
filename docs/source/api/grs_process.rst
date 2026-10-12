grs.grs_process
===============

.. automodule:: grs.grs_process
   :no-members:

:py:class:`~grs.grs_process.Process` runs the whole chain described in :doc:`../methods`:

* :py:meth:`~grs.grs_process.Process.execute` loads the Level-1 image (path to a ``.SAFE``, a Landsat
  Level-1 product or a NetCDF file, or an :py:class:`xarray.Dataset`), the CAMS data and the look-up
  tables, classifies the pixels and computes the Level-2 product, kept in memory in
  ``process.l2a`` (:py:class:`~grs.output.L2aProduct`);
* :py:meth:`~grs.grs_process.Process.write_output` writes it to NetCDF in ``odir``.

**Processing by blocks.** The image is cut into blocks of 512 x 512 pixels processed in parallel by a
:py:class:`multiprocessing.Pool` of ``ncpu`` processes (``processor.ncpu`` in ``config.yml``). The
results are written into shared arrays of the size of the full image: the memory needed grows with the
size of the image and the resolution (a Sentinel-2 tile at 10 m needs much more memory than at 20 or
60 m).

**Errors.** Reading errors of the input files (Level-1, CAMS, look-up tables, SURFWATER) are raised as
:py:class:`~grs.grs_exceptions.GRS_IO_Exception`.
