'''

Version history
==================

1.0.0:
    - grs for Sentinel2

1.1.0:
    - adaptation to Landsat 4, 5, 7, 8

1.1.1:
    - small changes

1.1.2:
    - output in Rrs unit

1.1.3
    - small changes

1.1.4:
    - set cloud mask; (compressed) netcdf4 output format

1.2.0:
    - load full data matrix from image instead of line by line pixel extraction (preparation for multipixel retrieval algorithm

1.2.1:
    - enable high latitude processing (change of DEM), option to process all pixels before masking "non-water" pixel

1.2.2:
    - Interpolation based on nearest-neighbor to keep tile-edge pixels.
    - Implementation of product.dispose to minimize memory usage in the jvm of snap

1.2.3:
    - compliant with version 8 of SNAP:
    - change output writing (now directly in NETCDF4, i.e., compressed, at the end of the process)
    - new utils get_subset

1.2.4:
    - option to load MAJA and WaterDetect mask and export masks in output file
    - option to process WaterDetect Water pixels only

1.3.0:
    - update for CAMS data (cds version):
    - new aod wavelengths and spectral ssa from 2018 onwards
    - adjustment for absorbing aerosol through ssa

1.3.1:
    - update CAMS data extraction from xarray and fix for longitude conventions

1.3.2:
    - Important improvement in LUT interpolation and access

1.3.3:
    - fix small bug for sunglint BRDF output

1.3.4:
    - add slope and shade from DEM

1.4.0:
    - process image by rectangular chunks

1.5.0:
    - change output parameter with addition of ndwi_nir and ndwi_swir

2.0.0:
    - remove snappy skeleton, simplification of the previous option, no more handling for Landsat

2.0.1:
    - fix for I/O

2.0.2:
    - memory optimization

2.0.3:
    - add surfwater

2.0.4:
    - adapt format for QGIS

2.0.5:
    - new output format compliant with SNAP "beam format"

2.1.0:
    - Major update on radiative transfer look-up tables and aerosol models

2.1.1:
    - package data

2.1.2:
    - fix for tiles straddling Greenwich meridian

2.1.3:
    - new cams automated loading and revised parameters

2.1.4:
    - change input/output feature for grs process
      (i.e., can be called as simple function within a script)
    - Change of output variable with addition of bitmask flags

2.1.5:
    - repackaging of some classes
    - include Landsat 8 and 9 imageries

2.1.6:
    - use GRSdriver instead s2driver
'''

__package__ = 'grs'
__version__ = '2.1.6'


from .acutils import Aerosol, Misc, Rasterization
from .cams import CamsProduct
from .auxdata import SensorData, AuxData
from .product import Product
from .output import L2aProduct
from .mask import Masking
from .grs_process import Process

import logging

#init logger
logger = logging.getLogger()

level = logging.getLevelName("INFO")
logger.setLevel(level)
