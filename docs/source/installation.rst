Installation
============

Look-up tables
--------------

The radiative transfer look-up tables (LUT) are not part of the package. Download them from
`grsdata <https://drive.google.com/drive/folders/1N0-FtW-PTPblR4z-82fFrUTekMd8e3Vz?usp=sharing>`__
and save them in a directory of your choice (``your_GRSDATA_PATH`` below).

Python environment
------------------

Python >= 3.9 is required (3.11 or 3.12 recommended). With conda:

.. code-block:: bash

   conda create -n grs python=3.12
   conda activate grs
   conda install -c conda-forge gdal eoreader cdsapi netCDF4 docopt xmltodict numba

Configuration
-------------

Set the paths of the LUTs in ``grs/config.yml`` before installing:

.. code-block:: yaml

   path:
     grsdata: 'your_GRSDATA_PATH'
     toa_lut: 'toa_lut_opac_wind_v3.nc'
     trans_lut: 'transmittance_lut_opac_wind_v3.nc'
     cams: 'your_CAMS_PATH'           # used when no CAMS file is given

   settings:
     aerosol_combination: [0, 0.5, 0, 0.5, 0., 0.]
     aerosol_models: ['ARCT_rh70', 'COAV_rh70', 'DESE_rh70', 'MACL_rh70', 'URBA_rh70', 'WASO_rh0']

   processor:
     ncpu: 8
     chunk: 512
     netcdf_engine: 'h5netcdf'

Install the package
-------------------

From the root of the repository:

.. code-block:: bash

   pip install .

Then check the command-line tool:

.. code-block:: bash

   grs -v

Build this documentation
------------------------

.. code-block:: bash

   pip install ".[docs]"
   sphinx-build -b html docs/source docs/build/html

The example notebooks are not executed: they are rendered with the outputs saved in them.
