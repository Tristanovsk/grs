from pathlib import Path
import yaml
import sys 

if __name__ == '__main__':


    if(len(sys.argv)>1):
        config_file=sys.argv[1]
    else:
        config_file="/app/grs/grs/config.yml"

    with open(config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)


    for key, value in data.items():
        if(value is not None and value!=''):
            config[key]=value 
            else:
            config[key]=None
            #return
            try:
                from grs import grs_process
                grs_process.process().execute(file_tbpd=config["file_tbpd"], outfile=config["outfile"], wkt=config["wkt"], 
                altitude=config["altitude"], aerosol=config["aerosol"],
                dem=config["dem"], aeronet_file=config["aeronet_file"], resolution=config["resolution"],
                aot550=config["aot550"], angstrom=config["angstrom"], unzip=config["unzip"], 
                untar=config["untar"], startrow=config["startrow"], angleonly=config["angleonly"])
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
