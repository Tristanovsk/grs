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
import geopandas as gpd
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

lev = 'L2GRS'

logdir = './tmp'

dirsat = '/datalake/'

l1cdir = {'s2': '/datalake/S2-L1C',
          'landsat': '/datalake/watcal/L8-L1-C2/'}
odir_root = {'s2': '/datalake/watcal/S2-L2GRS-aeronet/',
             'landsat': '/datalake/watcal/L8-L2GRS-C2-aeronet/'}

dem = True  # True
angleonly = False  # if true, grs is used to compute angle parameters only (no atmo correction is applied)
noclobber = True  # False #True
memory_safe = False  # True #
aeronet_file = 'no'
aerosol = 'cams'
# aerosol = 'user_model'
aot550 = 0.08  # used if aerosol = 'user_model'
angstrom = 1.6  # used if aerosol = 'user_model'
# Set radius of the subset-box in km
radius = 5
version='_v15'
# --------------------------------------------------------------------------------
basename = "tmp_grslist"
suffix = datetime.now().strftime("%y%m%d_%H%M%S")
tmp_file = "_".join([basename, suffix]) + '.tmp'


# -----------------
# function to get subset-box bounds
# -----------------
def get_buffer_box(p_lat, p_long, distance_km):
    # distance is d/2 of the square buffer around the point,
    # from center to corner;
    # find buffer width in meters
    buffer_width_m = (distance_km * 1000) / np.sqrt(2)

    # EPSG:4326 sets Coordinate Reference System to WGS84 to match input
    wgs84_pt_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([p_long], [p_lat], crs='EPSG:4326'))

    # find suitable projected coordinate system for distance
    utm_crs = wgs84_pt_gdf.estimate_utm_crs()
    # reproject to UTM -> create square buffer (cap_style = 3) around point -> reproject back to WGS84
    wgs84_buffer = wgs84_pt_gdf.to_crs(utm_crs).buffer(buffer_width_m, cap_style=3).to_crs('EPSG:4326')
    # wgs84_buffer.bounds returns bounding box as pandas dataframe,
    # .values[0] will extract first row as an array
    return wgs84_buffer.bounds


# --------------------------------------------------------------------------------

args_list = []
for idx, site in sites.iterrows():
    # load row of list file
    if site.iloc[0] == 0:
        continue
    name, start_date, end_date, sat, tile, resolution, flag, lat, lon = site.iloc[1:]
    if start_date == end_date:
        end_date = (datetime.strptime(end_date, '%Y-%m-%d').date() + timedelta(days=1)).__str__()
    sat = sat.lower()
    if sat == 'landsat':
        tile = '{:06d}'.format(tile)

    resolution = int(resolution)
    allpixels = True
    if flag == 1:
        allpixels = False

    # set subset (Area of interest)
    bounds = get_buffer_box(lat, lon, radius)
    # wkt = misc.wktbox(lon, lat,width=width,height=height)
    lon_min, lat_min, lon_max, lat_max = bounds.values[0]

    # ------------------
    # setting up loop on dates
    daterange = pd.date_range(start_date, end_date)
    for date in daterange:
        kwarg = ''

        # add subset option
        if True:
            kwarg += ' --longlat {:.4f},{:.4f},{:.4f},{:.4f}'.format(lon_max, lon_min, lat_max, lat_min)

        # ------------------
        # check if L1C exists and set directories / files
        subdir = opj(tile, '{:04d}'.format(date.year), '{:02d}'.format(date.month), '{:02d}'.format(date.day))
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
        else:
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
            suffix = '_cc' + cc_str + '_'+name + version
        else:
            suffix = name + version

        # ------------------
        # check / create output directory
        odir = opj(odir_root[sat.lower()], name)

        if not os.path.exists(odir):
            os.makedirs(odir)

        # force no suffix
        # name = ''
        # ------------------
        # get image basename for output naming
        basename = os.path.basename(l1c)
        outfile = misc.set_ofile(basename, odir=odir, suffix=suffix,level_name=lev)
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
            ancillary = 'cams_forecast'  # 'cams_reanalysis'

        # if 'cams' in aerosol:
        aerosol = ancillary

        command = 'grs ' + l1c + ' -o ' + outfile + ' --aerosol ' + aerosol + ' --dem --resolution ' + str(
            resolution) + kwarg

        # -------------------
        # fetch L2A MAJA
        # -------------------
        l2a_dir = opj(dirsat, 'S2-L2A-THEIA', subdir, '*')
        l2a_maja = glob.glob(opj(l2a_dir, 'S*MTD_ALL.xml'))
        if not l2a_maja:
            print(l2a_dir + ' not loaded on /datalake')
            # continue
            l2a_maja = None
        else:
            l2a_maja = l2a_maja[0]
            command += ' --maja ' + l2a_maja

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


command = pd.read_csv(tmp_file,sep='\t').values

with Pool(processes=ncore) as pool:
    pool.map(call, command)
    pool.close
