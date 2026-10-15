.. _api:

API reference
=============

The main classes are importable from the package itself:

.. code-block:: python

   import grs

   process = grs.Process()
   process.execute('S2A_MSIL1C_..._T31TFJ_....SAFE', cams_file='cams.nc', resolution=20,
                   odir='S2A_..._L2AGRS.zarr')
   process.write_output()

Processing chain
----------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.grs_process.Process

.. autosummary::
   :toctree: generated
   :nosignatures:

   grs.run.main

Input data
----------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.product.Product
   grs.cams.CamsProduct
   grs.auxdata.AuxData

Pixel classification
--------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.mask.Masking
   grs.mask.Settings

Atmospheric and sunglint correction
-----------------------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.grs_kernel.Kernel

Output product
--------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.output.L2aProduct

Configuration
-------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   grs.config.load_config
   grs.config.default_data_root
   grs.config.lut_path
   grs.config.check_luts
   grs.config.download_luts
   grs.config.init_user_config
   grs.config.main

Utilities
---------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.utils.GaseousTransmittance
   grs.utils.Gases
   grs.utils.Aerosol
   grs.utils.CamsParams
   grs.utils.Misc
   grs.utils.Rasterization
   grs.class_logger.ServiceLogger
   grs.grs_exceptions.GRSException
   grs.grs_exceptions.GRS_IO_Exception

Module overviews
----------------

Design notes of the main modules (dask processing, flag bits, output formats).

.. toctree::
   :maxdepth: 1

   api/grs_process
   api/grs_kernel
   api/mask
   api/output
