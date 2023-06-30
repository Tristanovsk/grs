''' Executable to process Sentinel-2 L1C images for aquatic environment

Usage:
  grs <input_file> [--cams_file file] [-o <ofile>] [--odir <odir>] [--resolution res] \
   [--levname <lev>] [--no_clobber] [--allpixels] [--surfwater file] [--snap_compliant]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input file to be processed

  --cams_file file     Absolute path of the CAMS file to be used (mandatory)

  -o ofile         Full (absolute or relative) path to output L2 image.
  --odir odir      Ouput directory [default: ./]
  --levname lev    Level naming used for output product [default: L2Agrs]
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --resolution=res  spatial resolution of the scene pixels
  --allpixels      force to process all pixels whatever they are masked (cloud, vegetation...) or not
  --surfwater file  Absolute path of the surfwater geotiff file to be used
  --snap_compliant  Export output to netcdf aligned with "beam" for ESA SNAP software

  Example:
      grs /data/satellite/S2/L1C/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE --cams_file /data/satellite/S2/cnes/CAMS/2022-07-31-cams-global-atmospheric-composition-forecasts.nc --resolution 60
  For CNES datalake:
      grs /datalake/S2-L1C/31TFJ/2022/07/31/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE --cams_file /datalake/watcal/ECMWF/CAMS/2022/07/31/2022-07-31-cams-global-atmospheric-composition-forecasts.nc --resolution 20

'''

import os, sys
from docopt import docopt
import logging
from . import __package__,__version__
from .grs_process import process


def main():
    args = docopt(__doc__, version=__package__ + '_' +__version__)
    print(args)

    file = args['<input_file>']
    lev = args['--levname']
    cams_file = args['--cams_file']
    surfwater_file = args['--surfwater']
    noclobber = args['--no_clobber']
    allpixels = args['--allpixels']
    resolution = int(args['--resolution'])
    snap_compliant = args['--snap_compliant']

    ##################################
    # File naming convention
    ##################################

    outfile = args['-o']
    if outfile == None:
        basename = os.path.basename(file)
        outfile = basename.replace('L1C', lev)
        outfile = outfile.replace('.SAFE', '').rstrip('/')

    odir = args['--odir']
    if odir == './':
        odir = os.getcwd()

    if not os.path.exists(odir):
        os.makedirs(odir)

    outfile = os.path.join(odir, outfile)

    if os.path.isfile(outfile) & noclobber:
        print('File ' + outfile + ' already processed; skip!')
        sys.exit()

    print('call grs_process for the following paramater. File:' +
                 file + ', output file:' + outfile +
                 ', cams_file:' + cams_file +
                 ', resolution:' + str(resolution))
    logging.info('call grs_process for the following paramater. File:' +
                 file + ', output file:' + outfile +
                 ', cams_file:' + cams_file +
                 ', resolution:' + str(resolution))
    process().execute(file, outfile, cams_file=cams_file, resolution=resolution,
                      allpixels=allpixels,surfwater_file=surfwater_file,snap_compliant=snap_compliant)
    return


if __name__ == "__main__":
    main()
