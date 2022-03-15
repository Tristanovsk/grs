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

    try:
        os.symlink(data['auxdata_path'], "/tmp/.snap/")
    except Exception as error:
        print(error)

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

    try:
        from grs import grs_process
        grs_process.process().execute(file=outfile, outfile=data["outfile"], wkt=wkt_rect, 
        altitude=data["altitude"], aerosol=data["aerosol"],
        dem=data["dem"], aeronet_file=data["aeronet_file"], resolution=data["resolution"],
        aot550=data["aot550"], angstrom=data["angstrom"], unzip=unzip, 
        untar=untar, startrow=data["startrow"])
    except Exception as inst:
        print('-------------------------------')
        print('error for file  ', inst, ' skip')
        print('-------------------------------')
        with open(config["logfile"], "a") as myfile:
            myfile.write('error during grs \n')

