''' Executable to process Sentinel-2 L1C images for aquatic environment

Usage:
  grs <input_file> [--cams_file file] [-o <ofile>] [--odir <odir>] [--resolution res] [--scale_aot factor]\
   [--opac_model name] [--levname <lev>] [--no_clobber] [--allpixels] [--surfwater file] [--dem_file file] [--snap_compliant]
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
  --dem_file file  Absolute path of the DEM geotiff file (already subset for the S2 tile)
  --scale_aot factor  scaling factor applied to CAMS aod550 raster
                    [default: 1]
  --opac_model name  Force the aerosol model (OPAC) to be 'name'
                     (choice: ['ARCT_rh70', 'COAV_rh70', 'DESE_rh70',
                     'MACL_rh70', 'URBA_rh70'])
  --snap_compliant  Export output to netcdf aligned with "beam" for ESA SNAP software

  Example:
      grs /data/satellite/S2/L1C/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE --cams_file /data/satellite/S2/cnes/CAMS/2022-07-31-cams-global-atmospheric-composition-forecasts.nc --resolution 60
  For CNES datalake:
      grs /work/datalake/S2-L1C/31TFJ/2023/06/16/S2B_MSIL1C_20230616T103629_N0509_R008_T31TFJ_20230616T111826.SAFE --cams_file /work/datalake/watcal/ECMWF/CAMS/2023/06/16/2023-06-16-cams-global-atmospheric-composition-forecasts.nc --odir /work/datalake/watcal/test --resolution 20 --dem_file /work/datalake/static_aux/MNT/COP-DEM_GLO-30-DGED_S2_tiles/COP-DEM_GLO-30-DGED_31TFJ.tif

'''

import os, sys
from docopt import docopt
import logging
from . import __package__, __version__
from .grs_process import Process


def main():
    args = docopt(__doc__, version=__package__ + '_' + __version__)
    print(args)

    file = args['<input_file>']

    lev = args['--levname']
    cams_file = args['--cams_file']
    surfwater_file = args['--surfwater']
    dem_file = args['--dem_file']
    noclobber = args['--no_clobber']
    allpixels = args['--allpixels']
    resolution = int(args['--resolution'])
    scale_aot = float(args['--scale_aot'])
    opac_model = args['--opac_model']
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

    if os.path.exists(outfile) & noclobber:
        print('File ' + outfile + ' already processed; skip!')
        sys.exit()

    logging.info('call grs_process for the following paramater. File:' +
                 file + ', output file:' + outfile +
                 f', cams_file:{cams_file}' +
                 ', resolution:' + str(resolution))


    process_ = Process()
    process_.execute(file,
                      ofile = outfile,
                      cams_file=cams_file,
                      resolution=resolution,
                      scale_aot=scale_aot,
                      opac_model=opac_model,
                      dem_file=dem_file,
                      allpixels=allpixels,
                      surfwater_file=surfwater_file,
                      snap_compliant=snap_compliant)
    process_.write_output()

    return


if __name__ == "__main__":
    main()
