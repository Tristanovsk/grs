import glob
import numpy as np
from osgeo import gdal

in_dir = '/Users/cole/Desktop/Sentinel/'

# Search directory for desired bands
red_file = glob.glob(in_dir + '**B04.jp2') # red band
nir_file = glob.glob(in_dir + '**B08.jp2') # nir band

# Define a function to calculate NDVI using band arrays for red, NIR bands
def ndvi(red, nir):
 return ((nir - red)/(nir + red))

# Open each band using gdal
red_link = gdal.Open(red_file[0])
nir_link = gdal.Open(nir_file[0])

# read in each band as array and convert to float for calculations
red = red_link.ReadAsArray().astype(np.float)
nir = nir_link.ReadAsArray().astype(np.float)

# Call the ndvi() function on red, NIR bands
ndvi2 = ndvi(red, nir)

# Create output filename based on input name
outfile_name = red_file[0].split('_B')[0] + '_NDVI.tif'

x_pixels = ndvi2.shape[0] # number of pixels in x
y_pixels = ndvi2.shape[1] # number of pixels in y

# Set up output GeoTIFF
driver = gdal.GetDriverByName('GTiff')

# Create driver using output filename, x and y pixels, # of bands, and datatype
ndvi_data = driver.Create(outfile_name,x_pixels, y_pixels, 1,gdal.GDT_Float32)

# Set NDVI array as the 1 output raster band
ndvi_data.GetRasterBand(1).WriteArray(ndvi2)

# Setting up the coordinate reference system of the output GeoTIFF
geotrans=red_link.GetGeoTransform() # Grab input GeoTranform information
proj=red_link.GetProjection() # Grab projection information from input file

# now set GeoTransform parameters and projection on the output file
ndvi_data.SetGeoTransform(geotrans)
ndvi_data.SetProjection(proj)
ndvi_data.FlushCache()
ndvi_data=None