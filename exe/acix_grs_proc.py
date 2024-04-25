'''
command to process images over the aeronet-oc sites
'''

import os, sys
import re
import numpy as np
import pandas as pd
import glob
import datetime
from multiprocessing import Pool
import googlemaps

# sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc, multi_process

misc = misc()

from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *

google_key='AIzaSyAZj_t2R7qQhkPYlzK1Et26DT4VA6zGGyE'
gmaps = googlemaps.Client(key=google_key)

# --------------------------------------------------------------------------------
# set parameters
aeronet = False
if aeronet:
    list_file = '/local/AIX/tristan.harmel/project/acix/AERONETOC_Matchups_List_harmel.xlsx'
    sites = pd.read_excel(list_file)
    list_file = '/local/AIX/tristan.harmel/project/acix/AERONETOC_Matchups_List_harmel.csv'
    sites = pd.read_csv(list_file)
else:
    list_file = '/local/AIX/tristan.harmel/project/acix/ACIXII_Aqua_PhaseII_scene_tile_IDs_harmel.xlsx'
    list_file = '/local/AIX/tristan.harmel/project/acix/ACIX_scene_tile_IDs_L1C_Updated_5_28_2019.xlsx'

    sites = pd.read_excel(list_file)

lev = 'L2grs'

logdir = './tmp'
idir_root = {'S2A': '/nfs/DD/S2/L1/ESA',
             'S2B': '/nfs/DD/S2/L1/ESA',
             'LANDSAT_5': '/nfs/DD/Landsat/L1/uncompressed',
             'LANDSAT_7': '/nfs/DD/Landsat/L1/uncompressed',
             'LANDSAT_8': '/nfs/DD/Landsat/L1/L1C1/untar' } #uncompressed'}

odir_root = {'S2A': '/nfs/DP/S2/L2/GRS/',
             'S2B': '/nfs/DP/S2/L2/GRS/',
             'LANDSAT_5': '/nfs/DP/Landsat/L2/GRS/',
             'LANDSAT_7': '/nfs/DP/Landsat/L2/GRS/',
             'LANDSAT_8': '/nfs/DP/Landsat/L2/GRS/'}
odir_sub = 'acix_test'
resolution = None
missions = ['all', 'S2', 'Landsat']
mission = missions[2]
if mission == 'Landsat':
    # number of images to process within one jpy virtual machine (i.e., for one load of snappy)
    Nimage = 4
    # number of processors to be used
    ncore = 17
else:
    Nimage = 2
    ncore = 6
    resolution = 10
download = False  # set to True if you want to download missing images
angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = False #True
memory_safe=False
aeronet_file = 'no'
aot550 = 0.1
angstrom = 0.5
full_tile = False
if full_tile:
    w, h = 200, 200
else:
    # set rectangle limits (width and height in meters) of the subscene to process
    w, h = 0.5, 0.5

fmissing = os.path.join(logdir, 'list_missing_files.txt')
fjunk = os.path.join(logdir, 'list_junk_files.txt')
with open(fmissing, 'w'):
    pass
# --------------------------------------------------------------------------------

args_list = []
sites_clean= sites.drop_duplicates(subset=sites.columns.values[1:])
#idx, site = list(sites.iterrows())[1]
for idx, site in sites_clean.iterrows():
    if site.iloc[0] != site.iloc[0]:
        continue

    if aeronet:# for acix phase I
        name, lat, lon, date_raw, time, basename = site.iloc[0:6]
        date = datetime.datetime.strptime(date_raw, '%d-%m-%Y') + datetime.timedelta(hours=time)
    else:# for ACIX phase II
        name, lat, lon, date_raw, time, basename, time_diff = site.iloc[0:7]
        time_diff = re.sub(':.*','',str(time_diff))
        try:
            date = datetime.datetime.strptime(date_raw, '%d-%m-%Y') + datetime.timedelta(hours=int(time_diff))
        except:
            date = date_raw + datetime.timedelta(hours=int(time_diff))

    # c = gmaps.elevation((lat, lon))
    # altitude = max(0,c[0].get('elevation')) #site.iloc[7]

    #############
    if date.year > 2016:
        aerosol = 'cams_forecast'
        # continue
    #############
    else:
        aerosol = 'cams_reanalysis'
    ancillary = aerosol
    sensor = misc.get_sensor(basename)
    if sensor == None:
        print('non standard image, not processed: ', basename)
        continue

    productimage = sensor[1]
    sat = sensor[2]

    # skip S2/Landsat if mission == Landsat/S2
    if (('Landsat' in productimage) & (mission == 'S2')) | (('S2' in productimage) & (mission == 'Landsat')):
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
        if 'S2' in productimage:
            # check if other version of the image exists
            sub = basename.split('_')
            new_basename = sub[0] + '_' + sub[1] + '_' + sub[2] + '*' + sub[5] + '*zip'
            f = glob.glob(os.path.join(idir, new_basename))
            if f != []:
                file = f[0]
                basename = os.path.basename(file)

        else: # for Landsat, need to untar
            # TODO untar L1C1/basename.tgz in uncompressed/basename
            tarfile = file.replace('uncompressed', 'L1C1') + '.tgz'

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

            # continue
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
    # if (os.path.isfile(outfile + ".dim") | os.path.isfile(outfile + ".nc")) & noclobber:
    if os.path.isfile(outfile + ".nc") & noclobber:
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
    print(site)
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
    if aeronet:
        altitude = site.iloc[7]
    else:
        # c = gmaps.elevation((lat, lon))
        #altitude = max(0, c[0].get('elevation'))
        altitude = 0

    args_list.append([file_tbp, outfile, wkt, altitude, aerosol, aeronet_file, ancillary, resolution, \
                      aot550, angstrom, memory_safe, unzip, untar, startrow, angleonly])

# reshape args_list to process several images (Nimage) on each processor
command = []
for args in misc.chunk(iter(args_list), Nimage):
    command.append([args, fjunk])

with Pool(processes=ncore) as pool:
    pool.map(multi_process().grs_call, command[0:10], 1)
    pool.close
