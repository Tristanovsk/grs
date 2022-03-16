from pathlib import Path
import yaml
import sys, os
import netCDF4 as nc
import geopandas as gpd
import logging
from logging.handlers import RotatingFileHandler

def shp2wkt(shapefile):
    print(shapefile)
    tmp = gpd.GeoDataFrame.from_file(shapefile)
    #tmp.to_crs(epsg=4326, inplace=True)
    return tmp.geometry.to_wkt().values[0]

if __name__ == '__main__':

    
    if(len(sys.argv)>1):
        config_file=sys.argv[1]
    else:
        config_file="/app/grs/global_config.yml"

    with open(config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # init logger
    logger = logging.getLogger()
    log_level="INFO"

    level = logging.getLevelName(log_level)
    logger.setLevel(level)

    # file handle
    log_file=data["logfile"]
    file_handler = RotatingFileHandler(log_file, 'a', 1000000, 1)
    formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d    %(levelname)s:%(filename)s::%(funcName)s:%(message)s',
                                  datefmt='%Y-%m-%dT%H:%M:%S')
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # stream handler
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    try:
        os.symlink(data['auxdata_path'], "/tmp/.snap/auxdata")
    except Exception as error:
        logger.error(error)

    os.environ['DATA_ROOT'] = data['data_root']
    os.environ['CAMS_PATH'] = data['cams_folder']

    for key, value in data.items():
        if(value is not None and value!=''):
            data[key]=value 
        else:
            data[key]=None
            

    if data["shapefile"]!=None:
        wkt = shp2wkt(shapefile)
    else:
        lonmin, lonmax = -180, 180
        latmin, latmax = -90,90#-21.13
        wkt = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"
    
    outfile=data['outfile']

    #if os.path.isfile(data['outfile'] + ".dim") & data['noclobber']:
    #    print('File ' + outfile + ' already processed; skip!')
    #    sys.exit(-1)

    #if os.path.isfile(outfile + ".dim.incomplete"):# & False:
    #    print('found incomplete File ' + outfile + '; skip!')
    
    file=data["input_file"]

    unzip = False
    if os.path.splitext(file)[-1] == '.zip':
        unzip = True
    
    untar = False
    if os.path.splitext(file)[-1] == '.tar':
        unzip = True 
        
    if data['outfile'] == None:

        basename=os.path.basename(file)
        outfile = basename.replace('L1C', "L2")
        outfile = outfile.replace('.SAFE', '').rstrip('/')
        outfile = outfile.replace('.zip', '').rstrip('/')
        outfile = outfile.replace('L1TP', "L2")
        outfile = outfile.replace('.txt', '').rstrip('/')

        outfile = os.path.join(data['output_dir'], outfile)
    else:
        outfile=data['outfile']

    dem=False
    if data["dem"]=="True":
        dem=True

    waterdetect_only=False
    if data["waterdetect_only"]=="True":
        waterdetect_only=True

    try:
        from grs import grs_process
        grs_process.process().execute(file=file, outfile=outfile, wkt=wkt, 
        altitude=data["altitude"], aerosol=data["aerosol"],
        dem=dem, aeronet_file=data["aeronet_file"],
        resolution=data["resolution"], aot550=data["aot550"], 
        angstrom=data["angstrom"], unzip=unzip, untar=untar, 
        startrow=data["startrow"], maja_xml=data["maja_xml"],
        waterdetect_file=data["waterdetect_file"], 
        waterdetect_only=waterdetect_only)
    except Exception as inst:
        logger.info('-------------------------------')
        logger.info('error for file  ', inst, ' skip')
        logger.info('-------------------------------')
        with open(data["logfile"], "a") as myfile:
            myfile.write('error during grs \n')

