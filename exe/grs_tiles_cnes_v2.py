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
from osgeo import gdal
from multiprocessing import Pool
import subprocess

# CNES lib for datalake managment
# import lxml
# from libamalthee import Amalthee

sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from exe.procutils import misc, multi_process

opj = os.path.join
misc = misc()

# --------------------------------------------------------------------------------
# set parameters
# data source to fill datalake
# amalthee = Amalthee('peps')
sitefile = 'exe/list/list_grs_cnes_obs2mod.csv'
sitefile = sys.argv[1]

# number of processors to be used
ncore = int(sys.argv[2])

# sitefile = 'exe/list/list_landsat_jegou.csv'
sites = pd.read_csv(sitefile)

dirsat = '/datalake/'
l1cdir = {'s2': '/datalake/S2-L1C',
          'landsat': '/datalake/watcal/L8-L1-C2/'}
odir_root = {'s2': '/datalake/watcal/S2-L2A-GRS/',
             'landsat': '/datalake/watcal/L8-L2GRS-C2/'}
cams_folder = '/datalake/watcal/ECMWF/CAMS'
noclobber = False # True  #True
lev = 'L2Agrs'
logdir = './tmp'

# --------------------------------------------------------------------------------
basename = "tmp_grslist"
suffix = datetime.now().strftime("%y%m%d_%H%M%S")
tmp_file = "_".join([basename, suffix]) + '.tmp'

args_list = []

for idx, site in sites.iterrows():
    # load row of list file
    if site.iloc[0] == 0:
        continue
    name, start_date, end_date, sat, tile, resolution, flag = site.iloc[1:]
    if start_date == end_date:
        end_date = (datetime.strptime(end_date, '%Y-%m-%d').date() + timedelta(days=1)).__str__()
    sat = sat.lower()
    if sat == 'landsat':
        tile = '{:06d}'.format(tile)

    resolution = int(resolution)
    allpixels = True
    if flag == 1:
        allpixels = False

    # ------------------
    # setting up loop on dates
    daterange = pd.date_range(start_date, end_date)
    for date in daterange:
        kwarg = ''

        # ------------------
        # check if L1C exists and set directories / files
        year = date.strftime('%Y')
        month = date.strftime('%m')
        day = date.strftime('%d')
        subdir = opj(tile, year, month, day)
        l1c_dir = opj(l1cdir[sat], subdir)
        if not os.path.exists(l1c_dir):
            continue

        if sat == 'landsat':
            l1c = glob.glob(opj(l1c_dir, 'L*.tar'))
            kwarg += ' --untar'
        else:
            l1c = glob.glob(opj(l1c_dir, 'S2*.SAFE'))
        if not l1c:
            # print(l1c_dir + ' not loaded on /datalake')
            continue

        l1c.sort()
        l1c = l1c[-1]

        # -------------------
        # get Metadata
        # -------------------
        # if bad definition of suffix
        if name != name:
            name = ''
        if sat == 's2':
            print('tt:' + l1c)
            filename = gdal.Open(opj(l1c, 'MTD_MSIL1C.xml'))
            metadata = filename.GetMetadata()
            cc = round(float(metadata['CLOUD_COVERAGE_ASSESSMENT']))
            cc_str = f'{cc:03}'
            print(metadata)
            suffix = '_cc' + cc_str + name
        else:
            suffix = name

        # ------------------
        # check / create output directory
        odir = opj(odir_root[sat.lower()], subdir)

        if not os.path.exists(odir):
            os.makedirs(odir)

        # force no suffix
        # name = ''
        # ------------------
        # get image basename for output naming
        basename = os.path.basename(l1c)
        outfile = misc.set_ofile(basename, odir=odir, suffix=suffix, level_name=lev)
        print(outfile)
        if noclobber & os.path.exists(outfile):
            print(basename, ' already processed; skipping image; set noclobber as False to force processing')
            continue

        sensor = misc.get_sensor(basename)
        if sensor == None:
            print('non standard image, not processed: ', basename)
            continue

        # ------------------
        #  CAMS data selection
        # TODO check data availability from CAMS and update the dates below

        cams_file = opj(cams_folder, year, month, day,
                        date.strftime('%Y-%m-%d') + '-cams-global-atmospheric-composition-forecasts.nc')
        if (not os.path.exists(cams_file)):
            cams_file = opj(cams_folder, year, date.strftime('%Y-%m') +
                            '_month_cams-global-atmospheric-composition-forecasts.nc')

        command = 'grs ' + l1c + ' -o ' + outfile + ' --cams_file ' + cams_file + ' --resolution ' + str(
            resolution) + kwarg

        # -------------------
        # fetch L2A MAJA
        # -------------------
        # l2a_dir = opj(dirsat, 'S2-L2A-THEIA', subdir, '*')
        # l2a_maja = glob.glob(opj(l2a_dir, 'S*MTD_ALL.xml'))
        # if not l2a_maja:
        #     print(l2a_dir + ' not loaded on /datalake')
        #     # continue
        #     l2a_maja = None
        # else:
        #     l2a_maja = l2a_maja[0]
        #     command+=' --maja '+l2a_maja

        if allpixels:
            command += ' --allpixels '

        with open(tmp_file, 'a') as file:
            file.write(command + '\n')


# with Pool(processes=ncore) as pool:
#     print(pool.map(multi_process().grs_cnes, args_list, 1))
#     pool.close()
#     pool.join()

def call(command):
    print(command)
    pipeline_out = subprocess.call(command, stderr=subprocess.STDOUT, shell=True)
    return


command = pd.read_csv(tmp_file).values
#
# with Pool(processes=ncore) as pool:
#     pool.map(call, command)
#     pool.close
