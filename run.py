''' Executable to process L1C images from Sentinel-2 and Landsat mission series

Usage:
  grs <input_file> <sensor> [-o <file>] [--shape <shp>] \
   [--longlat <longmax,longmin,latmax,latmin> ] \
   [--altitude=alt] [--aerosol=DB] [--aeronet=<afile>] \
   [--aot550=aot] [--angstrom=ang] \
   [--resolution=res] [--no_clobber]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input MTL.txt file to be processed
  <sensor>         sensor type: S2A, S2B, LANDSAT_5, LANDSAT_7, LANDSAT_8
  -o <file>       Full (absolute or relative) path to output L2 image.
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --shape shp      Process only data inside the given shape
  --longlat <longmax,longmin,latmax,latmin>
                   Restrict ROI to long max, long min, lat max, lat min in decimal degree
                   [default: 180, -180, 90, -90]
  --altitude=alt   altitude of the scene to be processed in meters
                   [default: 0]
  --aerosol=DB     aerosol data base to use within the processing
                   DB: cams_forecast, cams_reanalysis, aeronet, user_model
                   [default: cams_reanalysis]
  --aeronet=afile  if `--aerosol` set to 'aeronet', provide aeronet file to use
  --aot550=aot     if `--aerosol` set to 'user_model', provide aot550 value to be used
                   [default: 0.1]
  --angstrom=ang     if `--aerosol` set to 'user_model', provide aot550 value to be used
                   [default: 1]
  --resolution=res spatial resolution of the scene pixels

'''

import numpy as np
import geopandas as gpd
from docopt import docopt

from grs.config import *


def shp2wkt(shapefile):
    tmp = gpd.GeoDataFrame.from_file(shapefile)
    tmp.to_crs(epsg=4326, inplace=True)
    return tmp.geometry.values[0].to_wkt()


def main():
    args = docopt(__doc__, version=__package__ + ' ' + VERSION)
    print(args)

    file = args['<input_file>']
    sensor = args['<sensor>']
    shapefile = args['--shape']
    if (args['--shape'] == None):
        lonmax, lonmin, latmax, latmin = np.array(args['--longlat'].rsplit(','), np.float)

    noclobber = args['--no_clobber']

    altitude = float(args['--altitude'])
    resolution = args['--resolution']
    aerosol = args['--aerosol']
    aeronet_file = 'no'
    if aerosol == 'aeronet':
        aeronet_file = args['--aeronet_file']
    print(args)

    print(file, sensor, shapefile, altitude, aerosol, noclobber, aeronet_file, resolution)
    if shapefile != None:
        wkt = shp2wkt(shapefile)
    else:
        wkt = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
              + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
              + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"

    from grs.grs_process import process
    process().execute(file, sensor, wkt, altitude=altitude, aerosol=aerosol, noclobber=noclobber,
                      gdm=None, aeronet_file=aeronet_file, resolution=resolution, indband=None,
                      aot550=args['--aot550'], angstrom=args['--angstrom'])
    return

    # # # file = "/DATA/sensor/SENTINEL2//Villefranche/L1C/S2A_MSIL1C_20170206T102211_N0204_R065_T32TLP_20170206T102733.SAFE" # /DATA/sensor/SENTINEL2/lucinda/L1C/S2A_OPER_PRD_MSIL1C_PDMC_20160728T024349_R016_V20160728T003036_20160728T003036.SAFE/"
    file = '/DATA/Satellite/SENTINEL2/test/dimitri/L1C/S2A_MSIL1C_20170903T090551_N0205_R050_T35TLE_20170903T090723.SAFE/'
    file = '/DATA/Satellite/SENTINEL2/test/thierry/L1C/S2A_OPER_PRD_MSIL1C_PDMC_20160330T050601_R108_V20160326T103406_20160326T103406.SAFE/'
    file = '/home/harmel/Dropbox/rus/S2/gernez/L1C/S2A_MSIL1C_20170210T082051_N0204_R121_T33HYD_20170210T083752.SAFE/'
    sensor = 'S2A'
    lonmax = 180
    lonmin = -180
    latmax = 90
    latmin = -90
    altitude = 0
    file = '/DATA/Satellite/LANDSAT/L8/L1TP/LC08_L1TP_196029_20180914_20180928_01_T1/LC08_L1TP_196029_20180914_20180928_01_T1_MTL.txt'
    sensor = 'LANDSAT_8'

    file = '/DATA/Satellite/LANDSAT/L7/L1TP/LE07_L1TP_195029_20180915_20180915_01_RT/LE07_L1TP_195029_20180915_20180915_01_RT_MTL.txt'
    sensor = 'LANDSAT_7'

    file = '/DATA/Satellite/LANDSAT/L5/LT05_L1TP_191029_20111111_20161005_01_T1/LT05_L1TP_191029_20111111_20161005_01_T1_MTL.txt'
    sensor = 'LANDSAT_5'

    shapefile = '/DATA/OBS2CO/vrac/shape/SPO04/SPO04.shp'

    altitude = 780  # 0  # 234.
    aeronet_file = "/DATA/AERONET/V3/Thessaloniki_aod_V3.lev15"

    resolution = None
    aerosol = 'cams_forecast'
    noclobber = False
    if shapefile != None:
        wkt = shp2wkt(shapefile)
    else:
        wkt = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
              + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
              + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"

    from grs.grs_process import process
    process().execute(file, sensor, wkt, altitude=altitude, aerosol='cams_forecast', noclobber=noclobber,
                      gdm=None, aeronet_file=None, resolution=resolution, indband=None)


if __name__ == "__main__":
    main()
