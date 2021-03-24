'''
command to process images over the aeronet-oc sites

example:
python3 exe/grs_from_list.py exe/List_images_grs_template.csv
'''

import os, sys
import numpy as np
import pandas as pd
import glob
from datetime import datetime, timedelta
from multiprocessing import Pool

sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from exe.procutils import misc, multi_process

opj = os.path.join
misc = misc()

# --------------------------------------------------------------------------------
# set parameters

# number of processors to be used
ncore = 4

lev = 'L2grs'

logdir = './tmp'

dirsat = '/datalake/'

l1cdir = {'s2': '/datalake/S2-L1C',
          'landsat': '/datalake/L8-L1C/'}
odir_root = {'s2': '/datalake/watcal/S2-L2GRS/',
             'landsat': '/datalake/watcal/L8-L2GRS/'}

angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = True  # False #True
memory_safe = False  # True #
aeronet_file = 'no'
aerosol = 'cams'
# aerosol = 'user_model'
aot550 = 0.08  # used if aerosol = 'user_model'
angstrom = 1.6  # used if aerosol = 'user_model'
allpixels = False  # True

# sitefile = sys.argv[1]  # 'exe/List_images_grs_template.csv'
sitefile = 'exe/list_grs_cnes_template.csv'
sites = pd.read_csv(sitefile)
# --------------------------------------------------------------------------------

args_list = []

for idx, site in sites.iterrows():
    # load row of list file
    if site.iloc[0] == 0:
        continue
    name, start_date, end_date, sat, tile, resolution = site.iloc[1:]
    if start_date == end_date:
        end_date = (datetime.strptime(end_date, '%Y-%m-%d').date() + timedelta(days=1)).__str__()
    sat = sat.lower()
    resolution = int(resolution)

    # ------------------
    # setting up loop on dates
    daterange = pd.date_range(start_date, end_date)
    for date in daterange:
        # TODO add AMALTHEE call to load requested data

        # ------------------
        # check if L1C exists and set directories / files
        subdir = opj(tile, '{:04d}'.format(date.year), '{:02d}'.format(date.month), '{:02d}'.format(date.day))
        l1c_dir = opj(l1cdir[sat], subdir)
        if not os.path.exists(l1c_dir):
            continue

        # TODO modify for landsat images
        l1c = glob.glob(opj(l1c_dir, 'S2*.SAFE'))
        if not l1c:
            #print(l1c_dir + ' not loaded on /datalake')
            continue
        else:
            l1c = l1c[0]
        l2a_dir = opj(dirsat, 'S2-L2A-THEIA', subdir, '*')
        l2a_maja = glob.glob(opj(l2a_dir, 'S*MTD_ALL.xml'))
        if not l2a_maja:
            print(l2a_dir + ' not loaded on /datalake')
            continue
        else:
            l2a_maja = l2a_maja[0]

        # TODO clean up files on HPC-CNES for easy access
        # For the moment waterdetect set as None
        wd_dir = opj(dirsat, 'S2-L2A-THEIA', subdir)
        waterdetect = None  # glob.glob(opj(wd_dir, '*.tif'))[0]

        # ------------------
        # check / create output directory
        odir = opj(odir_root[sat.lower()], subdir)

        if not os.path.exists(odir):
            os.makedirs(odir)

        # if bad definition of suffix
        if name != name:
            name = ''
        # force no suffix
        name = ''
        # ------------------
        # get image basename for output naming
        basename = os.path.basename(l1c)
        outfile = misc.set_ofile(basename, odir=odir, suffix=name)
        print(outfile)

        sensor = misc.get_sensor(basename)
        if sensor == None:
            print('non standard image, not processed: ', basename)
            continue

        # ------------------
        #  CAMS data selection
        if date.year > 2016:
            ancillary = 'cams_forecast'
        else:
            ancillary = 'cams_reanalysis'

        if 'cams' in aerosol:
            aerosol = ancillary

        # skip if incomplete (enables multiprocess)
        if os.path.isfile(outfile + ".dim.incomplete"):  # & False:
            print('found incomplete File ' + outfile + '; skipped!')
            continue

        args_list.append([l1c, outfile, aerosol, aeronet_file, ancillary, resolution, \
                          l2a_maja, waterdetect, \
                          aot550, angstrom, memory_safe, allpixels, angleonly])
command = []
for args in misc.chunk(iter(args_list), 1):
    command.append(args)

with Pool(processes=ncore) as pool:
    pool.map(multi_process().grs_cnes, command,1)
    pool.close
