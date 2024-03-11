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
import issues.lut_interpol.read_lut as read_lut

#file='/sat_data/satellite/sentinel2/L1C/31TFJ/S2A_MSIL1C_20201004T104031_N0209_R008_T31TFJ_20201004T125253.SAFE'
file='/datalake/S2-L1C/21MYT/2022/01/14/S2A_MSIL1C_20220114T141051_N0301_R110_T21MYT_20220114T155152.SAFE'



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

#####################################
# LOAD LUT FOR ATMOSPHERIC CORRECTION
#####################################
# load lut
#l2h.lutfine='/home/harmel/DATA/projet/grs/grsdata/LUT/S2A/lut_osoaa_band_integrated_aot0.01_aero_rg0.10_sig0.46_ws2_pressure1015.2.nc'
#l2h.lutcoarse='/home/harmel/DATA/projet/grs/grsdata/LUT/S2A/lut_osoaa_band_integrated_aot0.01_aero_rg0.80_sig0.60_ws2_pressure1015.2.nc'

print('loading lut...', l2h.lutfine)
lutf = acutils.lut(l2h.band_names)
lutc = acutils.lut(l2h.band_names)
lutf.load_lut(l2h.lutfine, indband)
lutc.load_lut(l2h.lutcoarse, indband)



l2h.rot = l2h.SensorData.rot
# normalization of Cext to get spectral dependence of fine and coarse modes
nCext_f = lutf.Cext / lutf.Cext550
nCext_c = lutc.Cext / lutc.Cext550
l2h.set_outfile('test.nc')
l2h.wkt, lonmin, lonmax, latmin, latmax = _utils.get_extent(l2h.Product)
l2h.crs = str(l2h.Product.getBand(l2h.band_names[0]).getGeoCoding().getImageCRS())

        ######################################
        #      Create output l2 product
        #          'l2_product'
        ######################################

#l2h.create_product()
# l2h.load_data()
# w, h = l2h.width, l2h.height
#
#
# def remove_na(arr):
#     return arr[~np.isnan(arr)]
#
# sza_= remove_na(np.unique(l2h.sza.round(1)))
# vza_= remove_na(np.unique(l2h.vza.round(1)))
# azi_= remove_na(np.unique(l2h.razi.round(0)))

sza_ = np.array([30,40,55,56,60])
vza_= np.array([1,2,5,6.6,13])
azi_= np.array([90,100,180])



plt.figure()
lutc.refl.isel(sza=-1,azi=1,vza=1).plot(x='wl',hue='aot',marker='o')
plt.figure()
lutf.refl.isel(sza=-1,azi=1,vza=1).plot(x='aot',hue='wl',marker='o')

f_ = lutf.refl.interp(azi=azi_).interp(vza=vza_).interp(sza=sza_)
plt.figure()
f_.isel(sza=-1,azi=-1,vza=1).plot(x='wl',hue='aot')

rlut = f_
f_.dims
aot_ = np.array(f_.aot, dtype=float)#,order='F')
read_lut.read_lut(*f_.shape,aot_, sza_, azi_, vza_, rlut[::-1])
read_lut.read_lut(6, 11, 5, 3, 5,aot_, sza_, azi_, vza_, rlut)

##################################
#
######## END
#
##################################




# szalut=range(10,70)
# sza=18.5
# r_ = 100
# for i in range(len(szalut)):
#     r__ = abs(szalut[i] - sza)
#     if (r__ > r_) :
#         isza = i - 1
#         break
#     else:
#         isza=i
#         r_ = r__
#
# print(isza,sza,szalut[isza])
