from pathlib import Path
import yaml
import sys
import os 

if __name__ == '__main__':

    if(len(sys.argv)>1):
        config_file=sys.argv[1]
    else:
        config_file="/app/grs/grs/global_config.yml"

    with open(config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    sys.symlink(data['auxdata_path'], "/tmp/.snap/")
    os.environ['DATA_ROOT'] = data['data_root']
    os.environ['CAMS_PATH'] = data['cams_folder']

    for key, value in data.items():
        if(value is not None and value!=''):
            data[key]=value 
            else:
            data[key]=None
            
    lonmin, lonmax = -180, 180
    latmin, latmax = -90,90#-21.13
    wkt_rect = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"

    if os.path.isfile(outfile + ".dim") & data['noclobber']:
            print('File ' + outfile + ' already processed; skip!')
            continue
        # skip if incomplete (enables multiprocess)
    if os.path.isfile(outfile + ".dim.incomplete"):# & False:
        print('found incomplete File ' + outfile + '; skip!')
        continue

    try:
        from grs import grs_process
        grs_process.process().execute(file_tbpd=data["file_tbpd"], outfile=data["outfile"], wkt=wkt_rect, 
        altitude=data["altitude"], aerosol=data["aerosol"],
        dem=data["dem"], aeronet_file=data["aeronet_file"], resolution=data["resolution"],
        aot550=data["aot550"], angstrom=data["angstrom"], unzip=data["unzip"], 
        untar=data["untar"], startrow=data["startrow"], angleonly=data["angleonly"])
    except:
        print('-------------------------------')
        print('error for file  ', file_tbp, ' skip')
        print('-------------------------------')
        with open(config["logfile"], "a") as myfile:
            myfile.write(file_tbp + ' error during grs \n')
        continue

