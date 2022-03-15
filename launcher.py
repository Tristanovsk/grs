from pathlib import Path
import yaml
import sys
import os 

if __name__ == '__main__':

    if(len(sys.argv)>1):
        config_file=sys.argv[1]
    else:
        config_file="/app/grs/grs/config.yml"

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
            #return
            try:
                from grs import grs_process
                grs_process.process().execute(file_tbpd=data["file_tbpd"], outfile=data["outfile"], wkt=data["wkt"], 
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
        # here sys.exit instead of "return" to terminate and close snappy and free memory
        #sys.exit()
        return
