GRS processor
=============

GRS (Glint Removal for Sentinel-2-like sensors) performs the atmospheric correction of high spatial resolution
multispectral images (Sentinel-2 MSI, Landsat 8 and 9 OLI) over inland and coastal waters: it removes the
contribution of the atmosphere and of the sunglint and gives the remote-sensing reflectance ``Rrs`` of the
water, with the aerosol optical thickness, the sunglint reflectance and quality flags, in a Zarr (or NetCDF)
product.

- Source code: https://github.com/CNES/GRSprocessor
- Documentation: https://cnes.github.io/GRSprocessor
- Look-up tables (LUTs): https://doi.org/10.5281/zenodo.23394202

Installation
------------

Python 3.11 or newer is needed:

.. code-block:: bash

   pip install GRSprocessor

The package is ``grs`` (``import grs``) and it installs two commands, ``grs`` (processing) and
``grs_config`` (configuration and LUTs).

The cloud masks use OmniCloudMask, which installs PyTorch. On a computer without a GPU, the CPU build of
PyTorch is much smaller; install it first with
``pip install torch --index-url https://download.pytorch.org/whl/cpu``.

Look-up tables
--------------

The radiative transfer LUTs (about 1 GB) are not in the package: they are downloaded from
`Zenodo <https://doi.org/10.5281/zenodo.23394202>`__ (checksums verified) the first time an image is
processed, into the user data folder (``~/.local/share/grs/lut`` on Linux). To download them beforehand,
or to check the configuration in use:

.. code-block:: bash

   grs_config --download   # download the missing LUTs
   grs_config              # print the configuration and check the LUTs

To keep the LUTs in another folder (e.g. shared by several users), set it once with
``grs_config --init /path/to/lut`` or with the environment variable ``GRS_DATA_ROOT``. Without internet
access (e.g. on the nodes of a cluster), download the LUTs beforehand and switch the automatic download off
with ``GRS_NO_DOWNLOAD=1``.

Ancillary data
--------------

GRS needs the CAMS atmospheric composition (aerosols, gases, pressure, wind) of the date of the image: a
NetCDF file of the Copernicus Atmosphere Monitoring Service, given with ``--cams_file``. A digital elevation
model (``--dem_file``) is optional.

Usage
-----

.. code-block:: bash

   grs S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE \
       --cams_file cams_forecast_2022-07.nc --resolution 20 --odir ./L2A

``grs -h`` lists the options (resolution, resolution of the aerosol retrieval, maximum cloud cover,
aerosol model, DEM, water mask...). From Python:

.. code-block:: python

   import grs

   process = grs.Process()
   process.execute('S2B_MSIL1C_...SAFE', cams_file='cams_forecast_2022-07.nc', resolution=20)
   process.odir = 'S2B_MSI_L2AGRS_....zarr'
   process.write_output()

The processing is lazy and done block by block; ``processor: ncpu`` of the configuration (8 by default)
sets the number of threads, and therefore the memory used (about 10 GB for a Sentinel-2 tile at 20 m).
