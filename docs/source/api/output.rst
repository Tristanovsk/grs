grs.output
==========

.. automodule:: grs.output
   :no-members:

:py:meth:`L2aProduct.to_netcdf <grs.output.L2aProduct.to_netcdf>` writes, in the output directory
``<odir>`` (created if needed), two NetCDF4 files named after the directory ``<name>``:

``<name>.nc``
   full-resolution product: ``Rrs`` (wl, y, x), ``BRDFg``, ``aot550``, ``sza``, ``vza``, ``raa``,
   ``dem``, ``flags``, ``mask`` (y, x), with the coordinate ``central_wavelength`` (true central
   wavelength of each band) and the georeferencing (``spatial_ref``).
``<name>_anc.nc``
   ancillary data on the CAMS grid (dimensions ``xc``, ``yc``): the CAMS variables and the total
   gaseous transmittance (dimension ``wl_all``).

The variables are compressed (zlib, level 5) and stored as scaled ``int16`` (e.g. ``Rrs`` and ``BRDFg``:
scale 1e-5, offset 0.3; ``aot550``: scale 0.001), unpacked automatically by xarray. Existing files with
the same name are replaced.

With ``snap_compliant=True``, ``Rrs`` is split into one variable per band (``Rrs_443``, ``Rrs_490``...)
with the ``beam`` metadata profile, so that the product can be opened in ESA SNAP.
