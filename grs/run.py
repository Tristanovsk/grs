''' Executable to process L1C images from Sentinel-2 and Landsat mission series

Usage:
  grs <input_file> [--grs_a] [--sensor <sensor>] [-o <ofile>] [--odir <odir>] [--shape <shp>] [--wkt <wktfile>]\
   [--longlat <longmax,longmin,latmax,latmin> ] \
   [--altitude=alt] [--dem] [--aerosol=DB] [--aeronet <afile>] \
   [--aot550=aot] [--angstrom=ang] [--output param] [--resolution=res] \
   [--maja <maja_xml_file>] [--waterdetect <waterdetect_file>] [--waterdetect_pixels] \
   [--levname <lev>] [--no_clobber] [--memory_safe] [--unzip]\
   [--allpixels]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input file to be processed
  --grs_a          Apply the advanced (beta) version of GRS if set
  --sensor sensor Set the sensor type: S2A, S2B, LANDSAT_5, LANDSAT_7, LANDSAT_8
                    (by default sensor type is retrieved from input file name)
  -o ofile         Full (absolute or relative) path to output L2 image.
  --odir odir      Ouput directory [default: ./]
  --levname lev    Level naming used for output product [default: L2grs]
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --shape shp      Process only data inside the given shape
  --wkt wktfile    Process only data inside the given wkt file
  --longlat <longmax,longmin,latmax,latmin>
                   Restrict ROI to long max, long min, lat max, lat min in decimal degree
                   [default: 180, -180, 90, -90]
  --altitude=alt   altitude of the scene to be processed in meters
                   [default: 0]
  --dem            Use SRTM digital elevation model instead of --altitude (need internet connection)
  --aerosol=DB     aerosol data base to use within the processing
                   DB: cams_forecast, cams_reanalysis, aeronet, user_model
                   [default: cams_forecast]
  --aeronet=afile  if `--aerosol` set to 'aeronet', provide aeronet file to use
  --aot550=aot     if `--aerosol` set to 'user_model', provide aot550 value to be used
                   [default: 0.1]
  --angstrom=ang     if `--aerosol` set to 'user_model', provide aot550 value to be used
                   [default: 1]
  --maja maja_xml_file   use of mask from MAJA L2A images, path to xml ID of the L2A image
  --waterdetect waterdetect_file  use of water mask from waterdetect algorithm,
                    path to the appropriate WaterDetect data file
  --output param   set output unit: 'Rrs' or 'Lwn' [default: Rrs]
  --resolution=res  spatial resolution of the scene pixels
  --unzip          to process zipped images seamlessly
  --memory_safe    use generic resampler instead of S2resampler to save memory
                   (induces loss in angle resolution per pixel for S2)
  --allpixels      force to process all pixels whatever they are masked (cloud, vegetation...) or not
  --waterdetect_pixels  if waterdetect file is provided, process only the pixels masked as "water"
'''

import netCDF4 as nc # imported here to avoid conflicts on mistraou
import numpy as np
import geopandas as gpd
from docopt import docopt

from .config import *


def shp2wkt(shapefile):
    print(shapefile)
    tmp = gpd.GeoDataFrame.from_file(shapefile)
    #tmp.to_crs(epsg=4326, inplace=True)
    return tmp.geometry.values[0].to_wkt()


def main():
    args = docopt(__doc__, version=__package__ + ' ' + VERSION)
    print(args)

    file = args['<input_file>']
    grs_a = args['--grs_a']
    lev = args['--levname']
    if grs_a and lev == "L2grs":
        lev = "L2grsa"
    sensor = args['--sensor']
    shapefile = args['--shape']
    if (args['--shape'] == None):
        lonmax, lonmin, latmax, latmin = np.array(args['--longlat'].rsplit(','), np.float)

    noclobber = args['--no_clobber']
    unzip = args['--unzip']
    memory_safe = args['--memory_safe']
    altitude = float(args['--altitude'])
    dem = args['--dem']
    allpixels = args['--allpixels']
    resolution = args['--resolution']
    aerosol = args['--aerosol']
    aot550 = float(args['--aot550'])
    angstrom = float(args['--angstrom'])
    aeronet_file = None
    if aerosol == 'aeronet':
        aeronet_file = args['--aeronet']
    maja_xml = args['--maja']
    waterdetect_file = args['--waterdetect']
    waterdetect_only = args['--waterdetect_pixels']

    output=args['--output']

    ##################################
    # File naming convention
    ##################################

    outfile = args['-o']
    if outfile == None:

        basename=os.path.basename(file)
        outfile = basename.replace('L1C', lev)
        outfile = outfile.replace('.SAFE', '').rstrip('/')
        outfile = outfile.replace('.zip', '').rstrip('/')
        outfile = outfile.replace('L1TP', lev)
        outfile = outfile.replace('.txt', '').rstrip('/')
        # if 'S2' in sensor:
        #     outfile = basename.replace('L1C', lev)
        #     outfile = outfile.replace('.SAFE', '').rstrip('/')
        #     outfile = outfile.replace('.zip', '').rstrip('/')
        # elif 'LANDSAT' in sensor:
        #     outfile = basename.replace('L1TP', lev)
        #     outfile = outfile.replace('.txt', '').rstrip('/')
        # else:
        #     print('Not recognized sensor, please try again!')
        #     sys.exit()
    odir = args['--odir']
    if odir == './':
        odir = os.getcwd()
    outfile = os.path.join(odir, outfile)

    if os.path.isfile(outfile + ".nc") & noclobber:
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
    process().execute(file, outfile, wkt, grs_a= grs_a, sensor=sensor, altitude=altitude, aerosol=aerosol,
                      dem=dem, aeronet_file=aeronet_file, resolution=resolution,
                      maja_xml=maja_xml, waterdetect_file=waterdetect_file, waterdetect_only=waterdetect_only,
                      aot550=aot550, angstrom=angstrom, output=output, allpixels=allpixels, memory_safe=memory_safe,
                      unzip=unzip)
    return


if __name__ == "__main__":
    main()
