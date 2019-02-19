'''
command to process images over the aeronet-oc sites
'''

import os, sys
import numpy as np
import pandas as pd
import glob
import datetime



sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc

misc = misc()
from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *



#--------------------------------------------------------------------------------
# set parameters
sitefile = '/local/AIX/tristan.harmel/project/acix/AERONETOC_Matchups_List.xlsx'
lev = 'L2grs'
aerosol = 'cams_forecast'
odir_root = {'S2A':'/nfs/DP/S2/L2/GRS/',
             'S2B':'/nfs/DP/S2/L2/GRS/',
             'LANDSAT_5':'/nfs/DP/Landsat/L2/GRS/',
             'LANDSAT_7':'/nfs/DP/Landsat/L2/GRS/',
             'LANDSAT_8':'/nfs/DP/Landsat/L2/GRS/'}
odir_sub = 'acix'
Nimage = 1
noclobber = True
resolution = None
aeronet_file = 'no'
aot550 = 0.1
angstrom = 0.5
#--------------------------------------------------------------------------------


sites = pd.read_excel(sitefile)  # , sep=' ')

full_tile = True
if full_tile:
    w, h = 200, 200
else:
    w, h = 1, 1

lonmin, lonmax = -180, 180
latmin, latmax = -90, 90  # -21.13
wkt_rect = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
           + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
           + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"


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

    # get L1 image full path
    idir = dic[productimage]['path']
    file = os.path.join(idir, basename)
    #------------------------
    # input file naming
    # Warning: only .zip or .tgz are permitted
    file = file.replace('SAFE','zip')
    if productimage != 'S2_ESA':
        file = file.replace('.tgz','')
        file = file + '.tgz'
    print(sensor, name, file)

    #----------------------------------------------
    #  DOWNLOAD SECTION
    #----------------------------------------------
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
        download_image.mp_worker(command)

    #----------------------------------------------
    #  PROCESS SECTION
    #----------------------------------------------
    # check / create output directory
    odir = os.path.join(odir_root[sensor[0]], odir_sub)
    if not os.path.exists(odir):
        os.makedirs(odir)

    # load list of raw image files producing exception
    # during grs process (causes to be investigated)
    fjunk = os.path.join(odir, 'list_junk_files.txt')
    try:
        with open(fjunk) as f:
            junkfiles = f.read().splitlines()
    except:
        junkfiles = []
        with open(fjunk, 'w'):
            pass

    # set output file
    outfile = misc.set_ofile(basename, odir=odir, suffix='_'+name+'_GRS')

    # get area to be processed
    wkt = misc.wktbox(lon, lat, width=w, height=h)

    # skip file if listed in junkfile
    if file in junkfiles:
        print('File ' + outfile + ' listed in '+fjunk+'; skipped!'
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


    basename = os.path.basename(file)
    outfile, sensor = misc.set_ofile(basename, odir=odir)
    print('-------------------------------')
    print('call grs for ', outfile, sensor)
    print('-------------------------------')

    # check if already partially processed, if so get startrow value
    checksum = outfile + '.checksum'
    startrow = 0
    try:
        with open(checksum) as f:
            checkdata = f.read().splitlines()
        for s in checkdata:
            ss = s.split()
            if ss[0] == 'startrow':
                startrow = int(ss[1])
    except:
        pass

    try:
        from grs import grs_process
        grs_process.process().execute(file, outfile, wkt, altitude=altitude, aerosol=aerosol,
                                      gdm=None, aeronet_file=aeronet_file, resolution=resolution,
                                      aot550=aot550, angstrom=angstrom, unzip=unzip, untar=untar,
                                      startrow=startrow)
        #isuccess += 1
    except:
        # TODO note file name into log file
        print('-------------------------------')
        print('error for file  ', file, ' skip')
        print('-------------------------------')
        with open(fjunk, "a") as myfile:
            myfile.write(file + '\n')
        continue

######################################################
# TO BE CONTINUED

for idx, row in sites.iterrows():
    print(row)
    site = row.site

    imgs = glob.glob(idir + '*' + row.tile + '*')
    print(imgs.__len__())

    if False:
        wkt = misc.wktbox(row.lon, row.lat, width=w, height=h)
    else:
        wkt = wkt_rect

    # ----------------------
    # remove file names which won't be processed whatsoever:
    imgs_tbp = []
    for file in imgs:
        # TODO remove this part to process images from 2019
        if '201901' in file:
            continue

        if file in junkfiles:
            continue

        basename = os.path.basename(file)
        if 'incomplete' in basename:
            continue

        outfile, sensor = set_ofile(basename, odir=odir)
        print(outfile, sensor)
        # skip if already processed (the .dim exists)
        # if os.path.isfile(outfile + ".dim") & os.path.isdir(outfile + ".data") & noclobber:
        if os.path.isfile(outfile + ".dim") & noclobber:
            print('File ' + outfile + ' already processed; skip!')
            continue
        # skip if incomplete (enables multiprocess)
        if os.path.isfile(outfile + ".dim.incomplete"):  # & False:
            print('found incomplete File ' + outfile + '; skip!')
            continue

        imgs_tbp.append(file)
    if imgs_tbp == []:
        continue
    # ----------------------
    isuccess = 1
    for files in misc.chunk(iter(imgs_tbp), Nimage):

        for file in files:


        # comment next lines if you want to process the full series of images
        # warning: this can be very (prohibitively) consuming in memory,
        # recommended use: set Nimage parameter and call this code within different subprocess
        if isuccess > Nimage:
            sys.exit()
