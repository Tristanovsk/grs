.. _api:

API reference
=============

The main classes are importable from the package itself:

.. code-block:: python

   import grs

   process = grs.Process()
   process.execute('S2A_MSIL1C_..._T31TFJ_....SAFE', cams_file='cams.nc', resolution=20,
                   odir='S2A_MSIL2AGRS_..._T31TFJ_...')
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

Atmospheric correction utilities
--------------------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.acutils.GaseousTransmittance
   grs.acutils.Gases
   grs.acutils.Aerosol
   grs.acutils.CamsParams
   grs.acutils.Misc
   grs.acutils.Rasterization

Output product
--------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.output.L2aProduct

Logging and exceptions
----------------------

.. autosummary::
   :toctree: generated
   :template: autosummary/class.rst
   :nosignatures:

   grs.class_logger.ServiceLogger
   grs.grs_exceptions.GRSException
   grs.grs_exceptions.GRS_IO_Exception

Module overviews
----------------

Design notes of the main modules (processing by blocks, flag bits, output format).

.. toctree::
   :maxdepth: 1

   api/grs_process
   api/mask
   api/output
