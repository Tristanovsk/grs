'''
Version history:

1.0.0 grs for Sentinel2

1.1.0 adaptation to Landsat 4, 5, 7, 8

1.1.1

1.1.2 output in Rrs unit

1.1.3

1.1.4: set cloud mask; (compressed) netcdf4 output format

1.2.0: load full data matrix from image instead of line by line pixel extraction (preparation for multipixel retrieval algorithm

1.2.1: enable high latitude processing (change of DEM), option to process all pixels before masking "non-water" pixel

1.2.2:  - Interpolation based on nearest-neighbor to keep tile-edge pixels.
        - Implementation of product.dispose to minimize memory usage in the jvm of snap

1.2.3: compliant with version 8 of SNAP:
        - change output writing (now directly in NETCDF4, i.e., compressed, at the end of the process)
        - new utils get_subset

1.2.4: - option to load MAJA and WaterDetect mask and export masks in output file
       - option to process WaterDetect Water pixels only

1.3.0: update for CAMS data (cds version):
       - new aod wavelengths and spectral ssa from 2018 onwards
       - adjustment for absorbing aerosol through ssa

2.0.0: - process image by rectangular chunks
'''

__package__ = 'grs'
__version__ = '1.2.5'

from .config import *
from .acutils import aerosol, lut, misc, smac
from .auxdata import Aeronet, cams, sensordata
from .utils import info, utils
from .grs_process import process
