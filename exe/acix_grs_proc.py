'''
command to process images over the aeronet-oc sites
'''

import os, sys
import numpy as np
import pandas as pd
import glob
import datetime

from grs import grs_process
sys.path.extend([os.path.abspath(__file__)])
sys.path.extend(['/home/harmel/Dropbox/work/git/satellite_app/grs/exe'])
from procutils import misc
misc=misc()
import sid


sitefile = sys.argv[1]

odir = sys.argv[2]  # '/nfs/DP/S2/L2/GRS/aeronet-oc/netcdf'
if not os.path.exists(odir):
    os.makedirs(odir)


idir = '/nfs/DD/S2/L1/ESA/'
sitefile = '/DATA/Workshop/ACIX/AERONETOC_Matchups_List.xlsx'
sites = pd.read_excel(sitefile)  # , sep=' ')


full_tile = True
if full_tile:
    w, h = 200, 200
else:
    w, h = 1, 1

lev = 'L2grs'
aerosol = 'cams_forecast'
fjunk = os.path.join(odir, 'list_junk_files.txt')
Nimage = 1
noclobber = True
resolution = None
aeronet_file = 'no'
aot550 = 0.1
angstrom = 0.5

lonmin, lonmax = -180, 180
latmin, latmax = -90, 90  # -21.13
wkt_rect = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
           + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
           + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"

# load list of raw image files producing exception during grs process (causes to be investigated)
try:
    with open(fjunk) as f:
        junkfiles = f.read().splitlines()
except:
    junkfiles = []
    with open(fjunk, 'w'):
        pass

for idx, site in sites.iterrows():
    if site.iloc[0] != site.iloc[0]:
        continue
    name, lat, lon, date_raw, time, basename = site.iloc[0:6]

    # get date in pratical format
    date = datetime.datetime.strptime(date_raw,'%d-%m-%Y') + datetime.timedelta(hours=time)
    sensor = misc.get_sensor(basename)
    # check if image is already downloaded
    file = os.path.join(idir, basename)

    print(sensor, name, file)

    if not os.path.exists(file):
        sid.run.main(mode='console', productimage='S2_ESA', fromdate=date.strftime('%Y-%m-%d'), todate=date.strftime('%Y-%m-%d'),
                    detailedfile=None, fileBDDEMIL=None, shape=os.path.join(os.path.dirname(__file__), 'examples', 'ALL04.shp'),
  shapeID='CODE_LAC')


######################################################
# TO BE CONTINUED

for idx, row in sites.iterrows():
    print(row)
    site = row.site
    altitude = row.alt
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

            # unzip image file
            unzip = False
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
                grs_process.process().execute(file, outfile, wkt, altitude=altitude, aerosol=aerosol,
                                              gdm=None, aeronet_file=aeronet_file, resolution=resolution,
                                              aot550=aot550, angstrom=angstrom, unzip=unzip, startrow=startrow)
                isuccess += 1
            except:
                # TODO note file name into log file
                print('-------------------------------')
                print('error for file  ', file, ' skip')
                print('-------------------------------')
                with open(fjunk, "a") as myfile:
                    myfile.write(file + '\n')
                continue

        # comment next lines if you want to process the full series of images
        # warning: this can be very (prohibitively) consuming in memory,
        # recommended use: set Nimage parameter and call this code within different subprocess
        if isuccess > Nimage:
            sys.exit()
