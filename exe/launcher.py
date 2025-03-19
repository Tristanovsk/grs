from grs import class_logger
from grs import __version__
from grs.grs_process import Process
from pathlib import Path
import yaml
import sys
import logging
from datetime import datetime
from osgeo import gdal

from exe.procutils import misc

misc = misc()


def main():
    # read config and prepare environment
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = "/home/grs2/exe/global_config.yml"

    with open(config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # file handle
    logfile = Path(data['logfile'])
    log_folder = logfile.parent
    logfile.parent.mkdir(parents=True, exist_ok=True)
    class_logger.ServiceLogger(log_file=logfile, error_log=Path(log_folder, "error.log"), log_level=data['level'],
                               log_console=True)

    # get all config
    with open(data['hymotep_config'], 'r') as config_file:
        data.update(yaml.load(config_file, Loader=yaml.FullLoader))

    for key, value in data.items():
        if value is not None and value != '':
            data[key] = value
        else:
            data[key] = None
    file = Path(data["input_file"])

    if not file.exists():
        logging.error("Missing input file. Process stopped")
        exit(-1)
    if not data["cams_folder"]:
        logging.error("Missing CAMS folder. Process stopped")
        exit(-1)

    input_filename = file.name

    # Get CAMS file
    if Path(data['cams_folder']).is_file():
        cams_file = Path(data['cams_folder'])
    else:
        input_date = datetime.strptime(input_filename.split("_")[2], '%Y%m%dT%H%M%S').date()
        year = input_date.strftime('%Y')
        month = input_date.strftime('%m')
        day = input_date.strftime('%d')
        logging.info('Search for the daily CAMS file')
        cams_file = Path(data['cams_folder'], year, month, day,
                         input_date.strftime('%Y-%m-%d') + '-cams-global-atmospheric-composition-forecasts.nc')
        if not cams_file.exists():
            logging.info('No daily CAMS file found. Search for the monthly one')
            cams_file = Path(data['cams_folder'], year,
                             input_date.strftime('%Y-%m') + '_month_cams-global-atmospheric-composition-forecasts.nc')
    logging.info(f'CAMS file : {cams_file}')

    # Verify existence of inputs
    if not file.is_dir():
        logging.error("Input file doesn't exit. Process stopped")
        exit(-1)

    if not cams_file.is_file():
        logging.error("CAMS file doesn't exit. Process stopped")
        exit(-1)

    if data["surfwater_file"] and not Path(data["surfwater_file"]).is_file():
        logging.error("SurfWater file doesn't exit. Process stopped")
        exit(-1)

    # prepare outfile
    suffix = data["suffix"]
    if not suffix:
        suffix = f'_V{__version__}'

    output_dir = Path(data['output_dir'])

    outdir = misc.set_ofile(input_filename, odir=output_dir, level_name='L2AGRS', suffix=suffix)
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = Path(outdir,  f"{outdir.name}.nc")

    # skip if already processed
    if Path(outfile).is_file() & data["noclobber"]:
        logging.info(f'File {outfile} already processed; skip!')
        exit(-1)

    logging.info(f'call grs_process for the following parameters. '
                 f'File: {file} , '
                 f'output idirectory:  {outdir} , '
                 f'cams_file:  {cams_file} , '
                 f'surfwater_file: {data["surfwater_file"]} , '
                 f'resolution: {data["resolution"]} , '
                 f'allpixels: {data["allpixels"]} , '
                 f'snap_compliant: {data["snap_compliant"]})')

    # first check cloud cover (for S2, not implemented for Landsat)
    if 'MSIL1C' in input_filename:
        max_cc = data["max_cloud_cover"]
        f_ = gdal.Open(Path(file, 'MTD_MSIL1C.xml'))
        metadata = f_.GetMetadata()
        cc = float(metadata['CLOUD_COVERAGE_ASSESSMENT']) / 100
        if cc >= max_cc:
            logging.info('input file not processed since cloud cover {:.3f} is greater than {:.3f}'.format(cc, max_cc))
            return

    try:
        process_ = Process()
        process_.execute(file,
                         odir=outdir,
                         cams_file=cams_file,
                         resolution=data["resolution"],
                         scale_aot=data["scale_aot"],
                         opac_model=data["opac_model"],
                         dem_file=data["dem_file"],
                         allpixels=data["allpixels"],
                         surfwater_file=data["surfwater_file"],
                         snap_compliant=data["snap_compliant"])
        process_.write_output()

    except Exception as inst:
        logging.error('-------------------------------')
        message = 'error for file  ' + str(inst) + ' skip'
        logging.error(message)
        logging.error('-------------------------------')
        logging.error('error during grs', exc_info=True)

    finally:
        # Close logger and get stats
        class_logger.get_instance().close()


if __name__ == '__main__':
    main()
