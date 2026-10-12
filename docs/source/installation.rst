Installation
============

Look-up tables
--------------

The radiative transfer look-up tables (LUT) are not part of the package. Download them from
`grsdata <https://drive.google.com/drive/folders/1N0-FtW-PTPblR4z-82fFrUTekMd8e3Vz?usp=sharing>`__
and save them in a directory of your choice (``your_GRSDATA_PATH`` below).

Python environment
------------------

Python >= 3.12 is required. With conda:

.. code-block:: bash

   conda create -n grs python=3.12
   conda activate grs
   conda install -c conda-forge gdal

The other dependencies (``requirements.txt``, including
`GRSdriver <https://github.com/Tristanovsk/GRSdriver>`__) are installed with the package below.

Configuration
-------------

Set the paths of the LUTs in ``grs/config.yml`` before installing:

.. code-block:: yaml

   path:
     grsdata: 'your_GRSDATA_PATH'
     toa_lut: 'toa_lut_opac_wind_light_v2.nc'
     trans_lut: 'transmittance_lut_opac_wind_light_v2.nc'
     cams: 'your_CAMS_PATH'

   processor:
     ncpu: 8                          # number of processes of the multiprocessing pool
     chunk: 512

.. note::

   In this version, ``path.cams`` and ``processor.chunk`` are not used: the CAMS file must be given
   for each image (``--cams_file``) and the block size of the processing is 512 pixels.

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
