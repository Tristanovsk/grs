'''
command to process images over the aeronet-oc sites
'''

import os, sys
import numpy as np
import pandas as pd
import glob
from datetime import datetime, timedelta
from multiprocessing import Pool

# sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc, multi_process

misc = misc()

from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *

# --------------------------------------------------------------------------------
# set parameters
sitefile = sys.argv[1] #'exe/List_images_grs_template.csv'
#sitefile = 'exe/List_images_grs_brazil.csv'

# number of images to process within one jpy v- timedelta(days=1))irtual machine (i.e., for one load of snappy)
Nimage = 1
# number of processors to be used
ncore = 2


sites = pd.read_csv(sitefile)
lev = 'L2grs'

logdir = './tmp'
idir_root = {'s2': '/nfs/DD/S2/L1/ESA',
             'landsat': '/nfs/DD/landsat/L1/uncompressed'}

odir_root = {'s2': '/nfs/DP/S2/L2/GRS/',
             'landsat': '/nfs/DP/Landsat/L2/GRS/'}

download = True #False  # set to True if you want to download missing images
angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = True
memory_safe= True
aeronet_file = 'no'
aot550 = 0.1
angstrom = 0.5

fmissing = os.path.join(logdir, 'list_missing_files.txt')
fjunk = os.path.join(logdir, 'list_junk_files.txt')
with open(fmissing, 'w'):
    pass
# --------------------------------------------------------------------------------

args_list = []

#idx, site = list(sites.iterrows())[1]
for idx, site in sites.iterrows():

    if site.iloc[0] == 0:
        continue
    name, start, end, sat, tile, lat, lon, h, w, altitude, resolution = site.iloc[1:12]

    if download:

        cloudmax = str(60)

        if 'S2' in sat:
            productimage = 'S2_ESA'
            mission = ''
            script = dic[productimage]['script']
            write = dic[productimage]['path']
            auth = dic[productimage]['auth']
            # due to modification of the nomemclature at the two following date,
            # downloading is splitted into four cases:

            tiledate1=datetime.strptime('2016-12-01', '%Y-%m-%d').date()
            tiledate2=datetime.strptime('2017-04-01', '%Y-%m-%d').date()
            date1=datetime.strptime(start, '%Y-%m-%d').date()
            date2=datetime.strptime(end, '%Y-%m-%d').date()

            if (date1 < tiledate1) & (date2 < tiledate2):
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, start, (tiledate1- timedelta(days=1)).__str__(), productimage]
                download_image.mp_worker(command)
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, tiledate1.__str__(), end, productimage]
                download_image.mp_worker(command)
            elif (date1 <= tiledate1) & (date2 >= tiledate2):
                print(start, tiledate1.__str__())
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, start, (tiledate1- timedelta(days=1)).__str__(), productimage]
                download_image.mp_worker(command)
                print(tiledate1.__str__(), tiledate2.__str__())
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, tiledate1.__str__(), (tiledate2- timedelta(days=1)).__str__(), productimage]
                download_image.mp_worker(command)
                print(tiledate2.__str__(), end)
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, tiledate2.__str__(), end, productimage]
                download_image.mp_worker(command)
            elif (date1 > tiledate1) & (date2 > tiledate2):
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, start, (tiledate2- timedelta(days=1)).__str__(), productimage]
                download_image.mp_worker(command)
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, tiledate2.__str__(), end, productimage]
                download_image.mp_worker(command)
            else:
                command = [script, lat, lon, write, auth, tile, mission, cloudmax, start, end, productimage]
                download_image.mp_worker(command)
        else:
            productimage = 'Landsat_USGS'
            mission = 'LC8'
            script = dic[productimage]['script']
            write = dic[productimage]['path']
            auth = dic[productimage]['auth']
            command = [script, lat, lon, write, auth, tile, mission, cloudmax, start, end, productimage]
            download_image.mp_worker(command)
    continue
    if name != name:
        name=''
    odir_sub = tile
    sat=sat.lower()
    resolution = int(resolution)
    # get date in pratical format
    start = datetime.strptime(start, '%Y-%m-%d') #+ datetime.timedelta(hours=time)
    end = datetime.strptime(end, '%Y-%m-%d')

    files = glob.glob(os.path.join(idir_root[sat] ,'*' + tile + '*'))
    print(files.__len__())

    for file in files:
        #------------------
        # get date and images within the given, date range
        basename = os.path.basename(file)
        if 's2' in sat:
            date = basename.split('_')[2].split('T')[0]

        else:
            try:
                file = glob.glob(file + '/*MTL.txt')[0]
            except:
                with open(fjunk, "a") as myfile:
                    myfile.write(file + ' image is incomplete or missing \n')
                continue
            date = basename.split('_')[3]
        date = datetime.strptime(date,'%Y%m%d')
        if (date < start) | (date > end):
            continue

        #------------------
        #  CAMS data selection
        if date.year > 2016:
            aerosol = 'cams_forecast'
        else:
            aerosol = 'cams_reanalysis'

        basename = os.path.basename(file)
        sensor = misc.get_sensor(basename)
        if sensor == None:
            print('non standard image, not processed: ', basename)
            continue

        # ------------------------
        # input file naming
        # Warning: only .zip images are permitted (for landsat, please use uncompressed images
        file = file.replace('SAFE', 'zip')
        print(sensor, name, file)


        # ----------------------------------------------
        #  PROCESS SECTION
        # ----------------------------------------------
        # check / create output directory
        odir = os.path.join(odir_root[sat.lower()], odir_sub)
        if not os.path.exists(odir):
            os.makedirs(odir)

        # load list of raw image files producing exception
        # during grs process (causes to be investigated)
        try:
            with open(fjunk) as f:
                junkfiles = f.read().splitlines()
        except:
            junkfiles = []
            with open(fjunk, 'w'):
                pass

        # set output file
        outfile = misc.set_ofile(basename, odir=odir, suffix='_' + name + '_GRS')

        # get area to be processed
        wkt = misc.wktbox(lon, lat, width=w, height=h)

        # skip file if listed in junkfile
        if file in junkfiles:
            print('File ' + outfile + ' listed in ' + fjunk + '; skipped!')
            continue

        # skip if already processed (the .dim exists)
        if (os.path.isfile(outfile + ".dim") | os.path.isfile(outfile + ".nc")) & noclobber:
            print('File ' + outfile + ' already processed; skipped!')
            continue

        # skip if incomplete (enables multiprocess)
        if os.path.isfile(outfile + ".dim.incomplete"):  # & False:
            print('found incomplete File ' + outfile + '; skipped!')
            # continue

        # unzip image file
        unzip = False
        untar = False
        if os.path.splitext(file)[-1] == '.zip':
            unzip = True
        if os.path.splitext(file)[-1] == '.tgz':
            untar = True

        print('-------------------------------')
        print('call grs for ', outfile, sensor)
        print('-------------------------------')

        # check if already partially processed, if so get startrow value
        startrow = 0
        # TODO double check checksum files and location for optimization
        if True:
            checksum = outfile + '.checksum'
            try:
                with open(checksum) as f:
                    checkdata = f.read().splitlines()
                for s in checkdata:
                    ss = s.split()
                    if ss[0] == 'startrow':
                        startrow = int(ss[1])
            except:
                pass
        # break

        args_list.append([file, outfile, wkt, altitude, aerosol, aeronet_file, resolution, \
                          aot550, angstrom, memory_safe, unzip, untar, startrow, angleonly])

# reshape args_list to process several images (Nimage) on each processor
command = []
for args in misc.chunk(iter(args_list), Nimage):
    command.append([args, fjunk])

with Pool(processes=ncore) as pool:
    pool.map(multi_process().grs_call, command, 1)
    pool.close