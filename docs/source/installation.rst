Installation
============

GRS needs **Python 3.11 or newer**. It is a pure Python package, distributed on PyPI as
``GRSprocessor``: all its dependencies (including `GRSdriver <https://github.com/Tristanovsk/GRSdriver>`__,
which reads the L1C images) are installed by pip, and the GDAL Python bindings are not needed.

Python environment
------------------

Install GRS in a dedicated environment, with venv:

.. code-block:: bash

   python3.12 -m venv .venv
   source .venv/bin/activate          # Windows: .venv\Scripts\activate
   python -m pip install --upgrade pip

or with conda (Python only; the packages are then installed with pip):

.. code-block:: bash

   conda create -n grs python=3.12 pip
   conda activate grs

PyTorch (computers without GPU)
-------------------------------

The cloud masks use OmniCloudMask, which installs PyTorch. On a computer without a GPU, install the CPU
build of PyTorch first, much smaller than the default (CUDA) one:

.. code-block:: bash

   python -m pip install torch --index-url https://download.pytorch.org/whl/cpu

Install GRS
-----------

From PyPI (the package is ``grs``, ``import grs``, and it installs the commands ``grs`` and
``grs_config``):

.. code-block:: bash

   python -m pip install GRSprocessor

or from the source code (e.g. for the development version):

.. code-block:: bash

   git clone https://github.com/CNES/GRSprocessor.git
   cd GRSprocessor
   python -m pip install .
   # to develop (editable install with the test and build tools):
   python -m pip install -e ".[dev]"

Check the installation and the configuration in use:

.. code-block:: bash

   grs -v
   grs_config

To upgrade: ``python -m pip install --upgrade GRSprocessor``.

Look-up tables
--------------

The radiative transfer look-up tables (LUTs, about 1 GB of NetCDF files) are not in the package: they
are downloaded from `Zenodo <https://doi.org/10.5281/zenodo.23394202>`__ (checksums verified) the first
time an image is processed, into the default LUT folder (``~/.local/share/grs/lut`` on Linux). To
download them beforehand and check them:

.. code-block:: bash

   grs_config --download
   grs_config

- To keep the LUTs in another folder (e.g. a copy shared by several users), write it in your user
  configuration file with ``grs_config --init /path/to/lut``, or set the ``GRS_DATA_ROOT`` environment
  variable.
- Without internet access (e.g. on the nodes of a cluster), download the LUTs beforehand and switch the
  automatic download off with ``GRS_NO_DOWNLOAD=1`` (or ``zenodo: auto_download: false`` in the user
  configuration).

The configuration is read at run time: no need to reinstall the package after changing it. See
:doc:`configuration`.

CAMS atmospheric data
---------------------

GRS uses the CAMS global atmospheric composition forecasts (aerosols, ozone, water vapour, pressure,
wind) of the date of the image, from the `Copernicus Atmosphere Data Store
<https://ads.atmosphere.copernicus.eu>`__ (ADS). Give a CAMS NetCDF file to the processing, or let GRS
download it: create an ADS account, accept the licence of the dataset *CAMS global atmospheric
composition forecasts*, and write your API key in ``~/.cdsapirc``:

.. code-block:: text

   url: https://ads.atmosphere.copernicus.eu/api
   key: <your-personal-access-token>

See the :doc:`CAMS tutorial <examples/auxdata>` for the CAMS archive and the download.

Build this documentation
------------------------

.. code-block:: bash

   python -m pip install ".[docs]"
   sphinx-build -b html docs/source docs/build/html

The example notebooks are not executed: they are rendered with the outputs saved in them.
