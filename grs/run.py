''' Executable to process Sentinel-2 L1C images for aquatic environment

Usage:
  grs <input_file> [--cams_file file] [--odir <odir>] [--resolution res] [--max_cloud_cover max_cc] [--scale_aot factor]\
   [--opac_model name] [--levname <lev>] [--no_clobber] [--allpixels] [--surfwater file] [--dem_file file]\
   [--suffix suffix] [--snap_compliant]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input file to be processed

  --cams_file file     Absolute path of the CAMS file to be used (mandatory)

  -o <odir>, --odir <odir>  Full (absolute or relative) path to output L2 image.
  --levname lev    Level naming used for output product [default: L2AGRS]
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --resolution=res  spatial resolution of the scene pixels [default: 60]
  --max_cloud_cover max_cc  Skip process if image level 1 cloud cover is greater than max_cc
                            in decimal number [default: 1]
  --allpixels      force to process all pixels whatever they are masked (cloud, vegetation...) or not
  --surfwater file  Absolute path of the surfwater geotiff file to be used
  --dem_file file  Absolute path of the DEM geotiff file (already subset for the S2 tile)
  --scale_aot factor  scaling factor applied to CAMS aod550 raster
                      [default: 1]
  --opac_model name  Force the aerosol model (OPAC) to be 'name'
                     (choice: ['ARCT_rh70', 'COAV_rh70', 'DESE_rh70',
                     'MACL_rh70', 'URBA_rh70'])
  --suffix suffix  A suffix to append to the output dir name to personalize a run default value is "V<GRS version>"
  --snap_compliant  Export output to netcdf aligned with "beam" for ESA SNAP software


  Example:
      grs /data/satellite/S2/L1C/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE --cams_file /data/cams/world/cams_forecast_2022-07.nc --resolution 60
  For CNES datalake:
      grs /work/datalake/S2-L1C/31TFJ/2023/06/16/S2B_MSIL1C_20230616T103629_N0509_R008_T31TFJ_20230616T111826.SAFE --cams_file /work/datalake/watcal/ECMWF/CAMS/2023/06/16/2023-06-16-cams-global-atmospheric-composition-forecasts.nc --odir /work/datalake/watcal/test --resolution 20 --dem_file /work/datalake/static_aux/MNT/COP-DEM_GLO-30-DGED_S2_tiles/COP-DEM_GLO-30-DGED_31TFJ.tif

'''

import logging

from pathlib import Path
from docopt import docopt
from osgeo import gdal

from . import class_logger
from . import __package__, __version__
from .grs_process import Process
from exe.procutils import misc

misc = misc()


def main():
    args = docopt(__doc__, version=__package__ + '_' + __version__)
    print(args)

    file = Path(args['<input_file>'])

    lev = args['--levname']
    cams_file = args['--cams_file']
    surfwater_file = args['--surfwater']
    dem_file = args['--dem_file']
    noclobber = args['--no_clobber']
    allpixels = args['--allpixels']
    resolution = int(args['--resolution'])
    max_cc = float(args['--max_cloud_cover'])
    scale_aot = float(args['--scale_aot'])
    opac_model = args['--opac_model']
    snap_compliant = args['--snap_compliant']
    suffix = args['--suffix']

    ##################################
    # File naming convention
    ##################################
    basename = file.name
    # first check cloud cover (for S2, not implemented for Landsat)
    if ('MSIL1C' in basename) and ('SAFE' in basename):
        f_ = gdal.Open(Path(file, 'MTD_MSIL1C.xml'))
        metadata = f_.GetMetadata()
        cc = float(metadata['CLOUD_COVERAGE_ASSESSMENT']) / 100
        if cc >= max_cc:
            logging.info('input file not processed since cloud cover {:.3f} is greater than {:.3f}'.format(cc, max_cc))
            return

    odir = args['--odir']
    if odir == './':
        odir = Path.cwd()

    if not suffix:
        suffix = f'_V{__version__}'

    outdir = misc.set_ofile(basename, odir=odir, level_name=lev, suffix=suffix)
    outdir.mkdir(parents=True, exist_ok=True)

    class_logger.ServiceLogger(log_file=str(Path(outdir, 'log_file.log')), error_log="error.log", log_level='INFO',
                               log_console=False)

    outfile = Path(outdir, outdir.name + ".nc")

    # skip if already processed
    if outfile.is_file() & noclobber:
        logging.info(f'File {outfile} already processed; skip!')
        exit(-1)

    logging.info(f'call grs_process for the following paramater. '
                 f'File: {file}, '
                 f'output directory: {odir}, '
                 f'cams_file:{cams_file}, '
                 f'resolution: {resolution}')

    try:
        process_ = Process()
        process_.execute(file,
                         odir=outdir,
                         cams_file=cams_file,
                         resolution=resolution,
                         scale_aot=scale_aot,
                         opac_model=opac_model,
                         dem_file=dem_file,
                         allpixels=allpixels,
                         surfwater_file=surfwater_file,
                         snap_compliant=snap_compliant)
        process_.write_output()
    except Exception:
        logging.error("Fatal error in grs_process", exc_info=True)
    return


if __name__ == "__main__":
    main()
