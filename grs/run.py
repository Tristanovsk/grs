''' Executable to process L1C images from Sentinel-2 and Landsat mission series

Usage:
  grs <input_file> <sensor> [-o <ofile>] [--odir <odir>] [--shape <shp>] [--wkt <wktfile>]\
   [--longlat <longmax,longmin,latmax,latmin> ] \
   [--altitude=alt] [--aerosol=DB] [--aeronet=<afile>] \
   [--aot550=aot] [--angstrom=ang] \
   [--resolution=res] [--levname <lev>] [--no_clobber]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input MTL.txt file to be processed
  <sensor>         sensor type: S2A, S2B, LANDSAT_5, LANDSAT_7, LANDSAT_8
  -o ofile         Full (absolute or relative) path to output L2 image.
  --odir odir       Ouput directory [default: ./]
  --levname lev    Level naming used for output product [default: L2grs]
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --shape shp      Process only data inside the given shape
  --wkt wktfile    Process only data inside the given wkt file
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

import netCDF4 as nc # imported here to avoid conflicts on mistraou
import numpy as np
import geopandas as gpd
from docopt import docopt

from .config import *


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
    aot550 = float(args['--aot550'])
    angstrom = float(args['--angstrom'])
    aeronet_file = 'no'
    if aerosol == 'aeronet':
        aeronet_file = args['--aeronet_file']

    print(args)

    ##################################
    # File naming convention
    ##################################
    outfile = args['-o']
    if outfile == None:
        lev = args['--levname']
        basename=os.path.basename(file)
        if 'S2' in sensor:
            outfile = basename.replace('L1C', lev)
            outfile = outfile.replace('.SAFE', '').rstrip('/')
            outfile = outfile.replace('.zip', '').rstrip('/')
        elif 'LANDSAT' in sensor:
            outfile = basename.replace('L1TP', lev)
            outfile = outfile.replace('.txt', '').rstrip('/')
        else:
            print('Not recognized sensor, please try again!')
            sys.exit()

    outfile = os.path.join(args['--odir'], outfile)

    if os.path.isfile(outfile + ".dim") & os.path.isdir(outfile + ".data") & noclobber:
        print('File ' + outfile + ' already processed; skip!')
        sys.exit()

    print(file, sensor, outfile, shapefile, altitude, aerosol, noclobber, aeronet_file, resolution)
    if shapefile != None:
        wkt = shp2wkt(shapefile)
    elif args['--wkt'] != None:
        with open(args['--wkt'], 'r') as f:
            wkt = f.read()
    else:
        wkt = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
              + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
              + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"

    from .grs_process import process
    # TODO add **kargs for optional arg like ancillary (should be connected to aerosol for cams choice of forecast or reannalysis
    process().execute(file, outfile, sensor, wkt, altitude=altitude, aerosol=aerosol,
                      gdm=None, aeronet_file=aeronet_file, resolution=resolution,
                      aot550=aot550, angstrom=angstrom)
    return


if __name__ == "__main__":
    main()
