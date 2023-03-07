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
ncore =int(sys.argv[2])

#sitefile = 'exe/list/list_landsat_jegou.csv'
sites = pd.read_csv(sitefile)

lev = 'L2grs'

logdir = './tmp'

dirsat = '/datalake/'

l1cdir = {'s2': '/datalake/S2-L1C',
          'landsat': '/datalake/watcal/L8-L1-C2/'}
odir_root = {'s2': '/datalake/watcal/S2-L2GRS/',
             'landsat': '/datalake/watcal/L8-L2GRS-C2/'}

dem=True #True
angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = True # False #True
memory_safe = False  # True #
aeronet_file = 'no'
aerosol = 'cams'
# aerosol = 'user_model'
aot550 = 0.08  # used if aerosol = 'user_model'
angstrom = 1.6  # used if aerosol = 'user_model'
allpixels = False  # True


# --------------------------------------------------------------------------------

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
        tile ='{:06d}'.format(tile)
    resolution = int(resolution)
    allpixels = True
    if flag == 1:
        allpixels = False

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
        if sat =='landsat':
            l1c = glob.glob(opj(l1c_dir, 'L*.tar'))

        l1c = glob.glob(opj(l1c_dir, 'S2*.SAFE'))

        if not l1c:
            # print(l1c_dir + ' not loaded on /datalake')
            continue
        else:
            l1c = l1c[-1]

        #-------------------
        # get Metadata
        #-------------------
        filename = gdal.Open(opj(l1c,'MTD_MSIL1C.xml'))
        metadata = filename.GetMetadata()
        cc=round(float(metadata['CLOUD_COVERAGE_ASSESSMENT']))
        cc_str=f'{cc:03}'
        print(metadata)

        #-------------------
        # fetch L2A MAJA
        #-------------------
        l2a_dir = opj(dirsat, 'S2-L2A-THEIA', subdir, '*')
        l2a_maja = glob.glob(opj(l2a_dir, 'S*MTD_ALL.xml'))
        if not l2a_maja:
            print(l2a_dir + ' not loaded on /datalake')
            #continue
            l2a_maja = None
        else:
            l2a_maja = l2a_maja[0]

        # TODO clean up files on HPC-CNES for easy access
        # For the moment waterdetect set as None
        wd_dir = opj(dirsat, 'S2-L2A-THEIA', subdir)
        waterdetect = None  # glob.glob(opj(wd_dir, '*.tif'))[0]
        #for David processing
        if "guimaraes" in sitefile:
            wd_dir = opj('/work/datalake/watcal/GRS/wd_masks', subdir)

            waterdetect = glob.glob(opj(wd_dir, '*.tif'))
            if len(waterdetect)==0:
                continue
            else:
                waterdetect =waterdetect[0]

        # ------------------
        # check / create output directory
        odir = opj(odir_root[sat.lower()], subdir)

        if not os.path.exists(odir):
            os.makedirs(odir)

        # if bad definition of suffix
        if name != name:
            name = ''
        # force no suffix
        #name = ''
        # ------------------
        # get image basename for output naming
        basename = os.path.basename(l1c)
        outfile = misc.set_ofile(basename, odir=odir, suffix='_cc'+cc_str+name)
        print(outfile)
        if noclobber & os.path.exists(outfile + '.nc'):
            print(basename, ' already processed; skipping image; set noclobber as False to force processing')
            continue

        sensor = misc.get_sensor(basename)
        if sensor == None:
            print('non standard image, not processed: ', basename)
            continue

        # ------------------
        #  CAMS data selection
        # TODO check data availability from CAMS and update the dates below
        if date.year > 2018:
            ancillary = 'cds_forecast'
        elif (date.year == 2018) and (date.month > 6):
            ancillary = 'cds_forecast'
        else:
            ancillary = 'cams_forecast' #'cams_reanalysis'

        #if 'cams' in aerosol:
        aerosol = ancillary

        # skip if incomplete (enables multiprocess)
        # if os.path.isfile(outfile + ".dim.incomplete"):  # & False:
        #     print('found incomplete File ' + outfile + '; skipped!')
        #     continue
        print([l1c, outfile, aerosol, aeronet_file, ancillary, resolution, \
                          dem,l2a_maja, waterdetect, \
                          aot550, angstrom, memory_safe, allpixels, angleonly])
        args_list.append([l1c, outfile, aerosol, aeronet_file, ancillary, resolution, \
                          dem,l2a_maja, waterdetect, \
                          aot550, angstrom, memory_safe, allpixels, angleonly])
# command = []
# for args in misc.chunk(iter(args_list), 1):
#     command.append(args)


with Pool(processes=ncore) as pool:
    print(pool.map(multi_process().grs_cnes, args_list, 1))
    pool.close()
    pool.join()
