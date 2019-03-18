'''
command to process images over the aeronet-oc sites
'''

import os, sys
import numpy as np
import pandas as pd
import glob
import datetime
from multiprocessing import Pool

# sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc, multi_process

misc = misc()

from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *

from sentinelsat import SentinelAPI

api = SentinelAPI('harmel', 'kroumir1')

# --------------------------------------------------------------------------------
# set parameters
sitefile = 'exe/List_images_grs_template.csv'
odir_sub = 'gernez'
sites = pd.read_csv(sitefile)
lev = 'L2grs'

logdir = './tmp'
idir_root = {'S2A': '/nfs/DD/S2/L1/ESA',
             'S2B': '/nfs/DD/S2/L1/ESA',
             'LANDSAT_5': '/nfs/DD/Landsat/L1/uncompressed',
             'LANDSAT_7': '/nfs/DD/Landsat/L1/uncompressed',
             'LANDSAT_8': '/nfs/DD/Landsat/L1/uncompressed'}

odir_root = {'S2A': '/nfs/DP/S2/L2/GRS/',
             'S2B': '/nfs/DP/S2/L2/GRS/',
             'LANDSAT_5': '/nfs/DP/Landsat/L2/GRS/',
             'LANDSAT_7': '/nfs/DP/Landsat/L2/GRS/',
             'LANDSAT_8': '/nfs/DP/Landsat/L2/GRS/'}

resolution = None
missions = ['all', 'S2', 'Landsat']
mission = missions[1]
if mission == 'Landsat':
    # number of images to process within one jpy virtual machine (i.e., for one load of snappy)
    Nimage = 4
    # number of processors to be used
    ncore = 17
else:
    Nimage = 1
    ncore = 2

download = True #False  # set to True if you want to download missing images
angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = True

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
    if name != name:
        name=''
    odir_sub = tile

    # get date in pratical format
    start = datetime.datetime.strptime(start, '%Y-%m-%d') #+ datetime.timedelta(hours=time)
    end = datetime.datetime.strptime(end, '%Y-%m-%d')
    # TODO modify for landsat
    files = glob.glob(os.path.join(idir_root['S2A'] ,'*' + tile + '*'))
    print(files.__len__())
    imgs = []
    for file in files:
        date = file.split('_')[2].split('T')[0]
        date = datetime.datetime.strptime(date,'%Y%m%d')
        if (date>= start) & (date <= end):
            imgs.append(file)
    if imgs == []:
        continue

    wkt = misc.wktbox(lon, lat,width=w,height=h)

    for file in imgs:
        date = file.split('_')[2].split('T')[0]
        date = datetime.datetime.strptime(date,'%Y%m%d')
        #############
        # CAMS data selection
        if date.year > 2016:
            aerosol = 'cams_forecast'
        else:
            aerosol = 'cams_reanalysis'

        basename = os.path.basename(file)
        sensor = misc.get_sensor(basename)
        if sensor == None:
            print('non standard image, not processed: ', basename)
            continue

        productimage = sensor[1]
        sat = sensor[2]

        # skip S2/Landsat if mission == Landsat/S2
        if (('Landsat' in productimage) & (mission == 'S2')) | (('S2' in productimage) & (mission == 'Landsat')):
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
        odir = os.path.join(odir_root[sensor[0]], odir_sub)
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

        if "Landsat" in productimage:
            try:
                file_tbp = glob.glob(file + '/*MTL.txt')[0]
            except:
                with open(fjunk, "a") as myfile:
                    myfile.write(file + ' image is incomplete or missing \n')
                continue
        else:
            file_tbp = file

        print('-------------------------------')
        print('call grs for ', outfile, sensor)
        print('-------------------------------')

        # check if already partially processed, if so get startrow value
        startrow = 0
        # TODO double check checksum files and location for optimization
        if False:
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

        args_list.append([file_tbp, outfile, wkt, altitude, aerosol, aeronet_file, resolution, \
                          aot550, angstrom, unzip, untar, startrow, angleonly])

# reshape args_list to process several images (Nimage) on each processor
command = []
for args in misc.chunk(iter(args_list), Nimage):
    command.append([args, fjunk])

with Pool(processes=ncore) as pool:
    pool.map(multi_process().grs_call, command, 1)
    pool.close
