'''
command to download images from a datasheet list
'''

import os, sys
import re
import pandas as pd
import glob
import datetime
import multiprocessing

# sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc

misc = misc()
from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *


list_file = '/local/AIX/tristan.harmel/project/acix/AERONETOC_Matchups_List_harmel.xlsx'
list_file = '/local/AIX/tristan.harmel/project/acix/ACIXII_Aqua_PhaseII_scene_tile_IDs_harmel.xlsx'
list_file = '/local/AIX/tristan.harmel/project/acix/ACIX_scene_tile_IDs_L1C_Updated_5_28_2019.xlsx'
sites = pd.read_excel(list_file)  # , sep=' ')

missions = ['all', 'S2', 'Landsat']
mission = missions[1]

# number of processors to be used
ncore = 12
if mission == 'S2':
    ncore = 2

command = []
list = []
_file=''
for idx, site in sites.iterrows():
    if site.iloc[0] != site.iloc[0]:
        continue
    name, lat, lon, date_raw, time, basename, time_diff = site.iloc[0:7]
    time_diff = re.sub(':.*','',str(time_diff))
    print(name, lat, lon, date_raw, time, basename, time_diff)
    # get date in pratical format
    try:
        date = datetime.datetime.strptime(date_raw, '%d-%m-%Y') + datetime.timedelta(hours=int(time_diff))
    except:
        date = date_raw + datetime.timedelta(hours=int(time_diff))
    sensor = misc.get_sensor(basename)
    if sensor == None:
        print('non standard image, not processed: ', basename)
        continue
    productimage = sensor[1]
    sat = sensor[2]

    # skip S2/Landsat if mission == Landsat/S2
    if (('Landsat' in productimage) & (mission == 'S2')) | (('S2' in productimage) & (mission == 'Landsat')):
        continue

    idir = dic[productimage]['path']
    file = os.path.join(idir, basename)
    # ------------------------
    # input file naming
    # Warning: only .zip or .tgz are permitted
    file = file.replace('SAFE', 'zip')
    if productimage != 'S2_ESA':
        file = file.replace('.tgz', '')
        file = file + '.tgz'

    # check if image is already downloaded
    if (not os.path.exists(file)) & (file != _file):
        _file=file
        print(sensor, name, file)
        list.append(file)
        fromdate = date.strftime('%Y-%m-%d')
        todate = datetime.datetime.strftime(date + datetime.timedelta(days=2), '%Y-%m-%d')
        cloudmax = str(100)

        script = dic[productimage]['script']
        write = dic[productimage]['path']
        auth = dic[productimage]['auth']
        if productimage == 'S2_ESA':
            auth = os.path.abspath('/local/AIX/tristan.harmel/git/sat/sid/auxdata/apihub_th.txt')
        tile = misc.get_tile(basename)
        command.append([script, lat, lon, write, auth, tile, sat, cloudmax, fromdate, todate, productimage])
        print(fromdate,' to ',todate, basename)
        print(site)
        # download image
        # download_image.mp_worker(command)

p = multiprocessing.Pool(ncore)
p.map(download_image.mp_worker, command)
p.close()
