'''
Main program
'''

from pathlib import Path


import os, shutil
import zipfile
import tarfile
import glob
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')

from esasnappy import ProductData, ProductIO
from grs import info,utils, acutils, AuxData


#file='/sat_data/satellite/sentinel2/L1C/31TFJ/S2A_MSIL1C_20201004T104031_N0209_R008_T31TFJ_20201004T125253.SAFE'
file='/datalake/S2-L1C/21MYT/2022/01/14/S2A_MSIL1C_20220114T141051_N0301_R110_T21MYT_20220114T155152.SAFE'
file='/datalake/S2-L1C/57WWS/2016/10/10/S2A_MSIL1C_20161010T012652_N0204_R074_T57WWS_20161010T012650.SAFE'




##################################
# Get sensor auxiliary data
##################################

print('Get sensor auxiliary data')
_utils = utils()

sensor = _utils.get_sensor(file)
sensordata = AuxData.sensordata(sensor)
resolution = sensordata.resolution
indband = sensordata.indband


print("Reading...")
print(file)
product = ProductIO.readProduct(file)

##################################
# Generate l2h object
##################################

l2h = info(product, sensordata)
l2h.headerfile = file
#l2h.product = _utils.s2_resampler(l2h.product, resolution=resolution)
l2h.Product = _utils.resampler(l2h.Product, resolution=resolution)


##################################
# GET IMAGE AND RASTER PROPERTIES
##################################
l2h.get_bands(l2h.band_names)
l2h.print_info()
l2h.aux = AuxData.cams()
l2h.get_product_info()
l2h.wkt, lonmin, lonmax, latmin, latmax = _utils.get_extent(l2h.Product)

##################################
## ADD ELEVATION BAND
##################################
high_latitude = (latmax >= 60) | (latmin <= -60)
l2h.get_elevation(high_latitude)

plt.imshow(l2h.elevation)