from grs import class_logger
from grs.grs_process import Process
import yaml
import sys, os
import logging
from datetime import datetime
from osgeo import gdal

sys.path.extend([os.path.abspath(__file__)])
from procutils import misc

misc = misc()


def main():
    # read config and prepare environment
    if (len(sys.argv) > 1):
        config_file = sys.argv[1]
    else:
        config_file = "/home/grs2/exe/global_config.yml"

    with open(config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # file handle
    log_folder = os.path.dirname(data['logfile'])
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    class_logger.ServiceLogger(log_file=data['logfile'], output_dir=log_folder, log_level=data['level'],
                               log_console=True)

    # get all config
    with open(data['hymotep_config'], 'r') as config_file:
        data.update(yaml.load(config_file, Loader=yaml.FullLoader))

    for key, value in data.items():
        if (value is not None and value != ''):
            data[key] = value
        else:
            data[key] = None
    file = data["input_file"]

    if not file:
        logging.error("Missing input file. Process stopped")
        exit(-1)
    if not data["cams_folder"]:
        logging.error("Missing CAMS folder. Process stopped")
        exit(-1)

    input_filename = os.path.basename(file.rstrip('/'))

    # Get CAMS file
    if os.path.isfile(data['cams_folder']):
        cams_file = data['cams_folder']
    else:
        input_date = datetime.strptime(input_filename.split("_")[2], '%Y%m%dT%H%M%S').date()
        year = input_date.strftime('%Y')
        month = input_date.strftime('%m')
        day = input_date.strftime('%d')
        logging.info('Search for the daily CAMS file')
        cams_file = os.path.join(
            data['cams_folder'],
            year,
            month,
            day,
            input_date.strftime('%Y-%m-%d') + '-cams-global-atmospheric-composition-forecasts.nc'
        )
        if (not os.path.exists(cams_file)):
            logging.info('No daily CAMS file found. Search for the monthly one')
            cams_file = os.path.join(
                data['cams_folder'],
                year,
                input_date.strftime('%Y-%m') + '_month_cams-global-atmospheric-composition-forecasts.nc'
            )
    logging.info('CAMS file : ' + cams_file)

    # Verify existence of inputs
    if not os.path.isdir(file):
        logging.error("Input file doesn't exit. Process stopped")
        exit(-1)

    if not os.path.isfile(cams_file):
        logging.error("CAMS file doesn't exit. Process stopped")
        exit(-1)

    if data["surfwater_file"] and not os.path.isfile(data["surfwater_file"]):
        logging.error("SurfWater file doesn't exit. Process stopped")
        exit(-1)

    # prepare outfile
    suffix = '_V' + str(data["chain_version"])
    output_dir = data['output_dir']
    if not output_dir:
        output_dir = "./"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    outdir = misc.set_ofile(input_filename, odir=output_dir, level_name='L2AGRS', suffix=suffix)
    basename = os.path.basename(outdir)
    outfile = outdir + "/" + basename + ".nc"

    # skip if already processed
    if os.path.isfile(outfile) & data["noclobber"]:
        logging.info('File ' + outfile + ' already processed; skip!')
        exit(-1)

    logging.info('call grs_process for the following parameters. File: ' +
                 file + ', output file: ' + outfile +
                 ', cams_file: ' + cams_file +
                 ', surfwater_file: ' + str(data["surfwater_file"]) +
                 ', resolution: ' + str(data["resolution"]) +
                 ', allpixels: ' + str(data["allpixels"]) +
                 ', snap_compliant: ' + str(data["snap_compliant"]))

    # first check cloud cover (for S2, not implemented for Landsat)
    if 'MSIL1C' in input_filename:
        max_cc = data["max_cloud_cover"]
        f_ = gdal.Open(os.path.join(file, 'MTD_MSIL1C.xml'))
        metadata = f_.GetMetadata()
        cc = float(metadata['CLOUD_COVERAGE_ASSESSMENT']) / 100
        if cc >= max_cc:
            logging.info('input file not processed since cloud cover {:.3f} is greater than {:.3f}'.format(cc, max_cc))
            return

    try:
        process_ = Process()
        process_.execute(file,
                         ofile=outfile,
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
