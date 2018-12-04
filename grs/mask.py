
import numpy as np
from skimage.filters import threshold_otsu

# existing mask :
# THEIA L2
# LANDSAT (4,5,7,8)
# dictionnaires condition des masks
# InfoMask={'Sentinel_2A' :{'NameTHEIA':"SENTINEL2A",  'waterpixel':'MaskCDRarray == 1'},
#         'Sentinel_2B' :{'NameTHEIA':"SENTINEL2B", 'waterpixel':'MaskCDRarray == 1'},
#         'LANDSAT_8':{'waterpixel':'(MaskCDRarray == 324)|(MaskCDRarray == 388)|(MaskCDRarray == 836)|(MaskCDRarray == 900)'},
#         'LANDSAT_7':{'waterpixel':'(MaskCDRarray == 68)|(MaskCDRarray == 132)'},
#         'LANDSAT_5':{'waterpixel':'(MaskCDRarray == 68)|(MaskCDRarray == 132)'},
#         'LANDSAT_4':{'waterpixel':'(MaskCDRarray == 68)|(MaskCDRarray == 132)'}}
# creating mask
# Otsu
# OBS2CO


class mask:
    def __init__(self):
        pass

    def water_detection(ndwi):
        """
        Very simple water detector based on Otsu thresholding method of NDWI.
        """
        otsu_thr = 1.0
        if len(np.unique(ndwi)) > 1:
            otsu_thr = threshold_otsu(ndwi)

        return ndwi>otsu_thr

