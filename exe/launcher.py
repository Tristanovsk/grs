from pathlib import Path
import yaml
import sys, os
import netCDF4 as nc
import geopandas as gpd
import logging
from logging.handlers import RotatingFileHandler

sys.path.extend([os.path.abspath(__file__)])
from procutils import misc
misc=misc()

def rename_file(file, outfile, outdir):

    if outfile == None:
        basename=os.path.basename(file)
        outfile = basename.replace('L1C', "L2GRS")
        outfile = outfile.replace('.SAFE', '').rstrip('/')
        outfile = outfile.replace('.zip', '').rstrip('/')
        outfile = outfile.replace('L1TP', "L2GRS")
        outfile = outfile.replace('.txt', '').rstrip('/')

    if outdir is None:
        return outfile
    else:
        return os.path.join(outdir, outfile)

def shp2wkt(shapefile):
    logging.info(f'{shapefile} is used')
    tmp = gpd.GeoDataFrame.from_file(shapefile)
    return tmp.geometry.to_wkt().values[0]

if __name__ == '__main__':

    #read config and prepare environment
    if(len(sys.argv)>1):
        config_file=sys.argv[1]
    else:
        config_file="/app/grs/exe/global_config.yml"

    with open(config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    try:
        os.symlink(data['auxdata_path'], "/tmp/.snap/auxdata")
    except Exception as error:
        logging.debug(error)

    os.environ['DATA_ROOT'] = data['data_root']
    os.environ['CAMS_PATH'] = data['cams_folder']


    from grs import grs_process

    for key, value in data.items():
        if(value is not None and value!=''):
            data[key]=value 
        else:
            data[key]=None
    file = data["input_file"] 
  
    if data["shapefile"] != None:
        wkt = shp2wkt(data["shapefile"])
    else:
        lonmin, lonmax = -180, 180
        latmin, latmax = -90,90#-21.13
        wkt = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"
    
    unzip = False
    if os.path.splitext(file)[-1] == '.zip':
        unzip = True
    
    untar = False
    if os.path.splitext(file)[-1] == '.tar':
        unzip = True 
        
    #outfile = rename_file(file, data["outfile"], data["output_dir"])

    dem=False
    if data["dem"]:
        dem=True

    waterdetect_only=False
    if data["waterdetect_only"]:
        waterdetect_only=True

    basename = os.path.basename(file)
    if 'incomplete' in basename:
        exite(-1)

    suffix='_'+str(data["chain_version"])+"_"+str(data["product_counter"])
    
    outfile = misc.set_ofile(file, odir=data['output_dir'], level_name='l2grs', suffix=suffix)
    logging.info(outfile)
    
    # skip if already processed (the .dim exists)
    if os.path.isfile(outfile + ".dim") & data["noclobber"]:
        logging.info('File ' + outfile + ' already processed; skip!')
        exit(-1)
    # skip if incomplete (enables multiprocess)
    if os.path.isfile(outfile + ".incomplete"):
        logging.info('found incomplete File ' + outfile + '; skip!')
        exit(-1)

    if os.path.isfile(outfile+".nc") & data["noclobber"]:
        logging.info('File ' + outfile + ' already processed; skip!')
        exit(-1)
   

    checksum = outfile+'.checksum'
    startrow=0
    try:
      with open(checksum) as f:
          checkdata = f.read().splitlines()
      for s in checkdata:
          ss=s.split()
          if ss[0]=='startrow':
              startrow=int(ss[1])
    except:
       pass


    try:
        grs_process.process().execute(file=file, outfile=outfile, wkt=wkt, 
        altitude=data["altitude"], aerosol=data["aerosol"],
        dem=dem, aeronet_file=data["aeronet_file"],
        resolution=data["resolution"], aot550=data["aot550"], 
        angstrom=data["angstrom"], unzip=unzip, untar=untar, 
        startrow=startrow, maja_xml=data["maja_xml"],
        waterdetect_file=data["waterdetect_file"], 
        waterdetect_only=waterdetect_only, memory_safe=data["memory_safe"], 
        angleonly=data["angleonly"], grs_a=data["grs_a"], output=data["output"])
    except Exception as inst:
        logging.info('-------------------------------')
        logging.info('error for file  ', inst, ' skip')
        logging.info('-------------------------------')
        with open(data["logfile"], "a") as myfile:
            myfile.write('error during grs \n')

