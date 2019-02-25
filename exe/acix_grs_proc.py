'''
command to process images over the aeronet-oc sites
'''

import os, sys
import numpy as np
import pandas as pd
import glob
import datetime
import multiprocessing

# sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc, multi_process

misc = misc()
from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *

# --------------------------------------------------------------------------------
# set parameters
sitefile = '/local/AIX/tristan.harmel/project/acix/AERONETOC_Matchups_List_harmel.xlsx'
lev = 'L2grs'
aerosol = 'cams_forecast'
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
odir_sub = 'acix'
missions=['all','S2','Landsat']
mission=missions[2]
# number of images to process within one jpy virtual machine (i.e., for one load of snappy)
Nimage = 2
# number of processors to be used
ncore = 5

download = False  # set to True if you want to download missing images
angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = True
resolution = None
aeronet_file = 'no'
aot550 = 0.1
angstrom = 0.5
full_tile = False
if full_tile:
    w, h = 200, 200
else:
    # set rectangle limits (width and height in meters) of the subscene to process
    w, h = 1, 1

fmissing = os.path.join(logdir, 'list_missing_files.txt')
fjunk = os.path.join(logdir, 'list_junk_files.txt')
with open(fmissing, 'w'):
            pass
# --------------------------------------------------------------------------------

args_list = []
sites = pd.read_excel(sitefile)
idx, site = list(sites.iterrows())[1]
for idx, site in sites.iterrows():
    if site.iloc[0] != site.iloc[0]:
        continue
    name, lat, lon, date_raw, time, basename = site.iloc[0:6]
    altitude = 0  # site.iloc[6]
    # get date in pratical format
    date = datetime.datetime.strptime(date_raw, '%d-%m-%Y') + datetime.timedelta(hours=time)

    sensor = misc.get_sensor(basename)
    if sensor == None:
        print('non standard image, not processed: ', basename)
        continue

    productimage = sensor[1]
    sat = sensor[2]

    # skip S2/Landsat if mission == Landsat/S2
    if (('Landsat' in productimage) & (mission=='S2')) | (('S2' in productimage) & (mission=='Landsat')):
        continue

    # get L1 image full path
    idir = idir_root[sensor[0]]
    file = os.path.join(idir, basename)
    # ------------------------
    # input file naming
    # Warning: only .zip images are permitted (for landsat, please use uncompressed images
    file = file.replace('SAFE', 'zip')
    print(sensor, name, file)

    # ----------------------------------------------
    #  DOWNLOAD SECTION
    # ----------------------------------------------
    # check if image is already downloaded
    if not os.path.exists(file):
        fromdate = date.strftime('%Y-%m-%d')
        todate = datetime.datetime.strftime(date + datetime.timedelta(days=1), '%Y-%m-%d')
        cloudmax = str(100)

        script = dic[productimage]['script']
        write = dic[productimage]['path']
        auth = dic[productimage]['auth']
        tile = misc.get_tile(basename)
        command = [script, lat, lon, write, auth, tile, sat, cloudmax, fromdate, todate, productimage]
        print(command)
        print('downloading image, please wait...')
        # download image
        if download:
            download_image.mp_worker(command)
        else:
            with open(fmissing, "a") as myfile:
                myfile.write(file + '\n')

            continue
    else:
        pass  # break
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
        continue

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
command=[]
for args in misc.chunk(iter(args_list), Nimage):
    command.append([args,fjunk])

p = multiprocessing.Pool(ncore)
p.map(multi_process().grs_call, command)


