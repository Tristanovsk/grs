from osgeo import gdal
import numpy as np
import rasterio
import xarray as xr
import richdem as rd

import matplotlib.pyplot as plt
import matplotlib

file='./test/results/S2B_MSIL2grs_20180927T103019_N0206_R108_T31TGK_20180927T143835.nc'
img = xr.open_dataset(file)

sza=np.radians(65)
sun_azi=np.radians(120)

dem = np.array(img.elevation)
dem= rd.rdarray(dem,no_data=np.nan)
slope= np.radians(rd.TerrainAttribute(dem, attrib="slope_degrees"))
aspect = rd.TerrainAttribute(dem, attrib="aspect")
aspect_rd = np.radians(aspect)
shaded = np.cos(sza) * np.cos(slope) + np.sin(sza)\
         * np.sin(slope) * np.cos(sun_azi - aspect_rd)

fig,axs=plt.subplots(2,2,figsize=(10,8))
axs=axs.ravel()
axs[0].imshow(dem,cmap=plt.cm.gist_earth)
axs[0].imshow(shaded,cmap=plt.cm.gray,alpha=0.5)
axs[1].imshow(slope,cmap=plt.cm.gray)
axs[2].imshow(aspect,cmap=plt.cm.Spectral_r)
axs[3].imshow(255*(shaded+1)/2,cmap=plt.cm.gray)

