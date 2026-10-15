Installation
============

Python environment
------------------

Python >= 3.9 is required (3.11 or 3.12 recommended). With conda:

.. code-block:: bash

   conda create -n grs python=3.12
   conda activate grs
   conda install -c conda-forge gdal eoreader cdsapi netCDF4 docopt xmltodict numba

Install the package
-------------------

From the root of the repository:

.. code-block:: bash

   pip install .

Then check the command-line tool:

.. code-block:: bash

   grs -v

Look-up tables
--------------

The radiative transfer look-up tables (LUTs) are not part of the package. Download them from Zenodo
into the default LUT folder (``~/.local/share/grs/lut`` on Linux) and check the configuration:

.. code-block:: bash

   grs_config --download
   grs_config

To keep the LUTs in another folder, first write it in your user configuration file with
``grs_config --init /path/to/lut``. See :doc:`configuration` for the other settings.

Build this documentation
------------------------

.. code-block:: bash

   pip install ".[docs]"
   sphinx-build -b html docs/source docs/build/html

The example notebooks are not executed: they are rendered with the outputs saved in them.
