'''
command to process images over the aeronet-oc sites
'''

import os, sys
import pandas as pd
import glob
import datetime
import multiprocessing

#sys.path.extend([os.path.abspath(__file__)])
sys.path.extend([os.path.abspath('exe')])
from procutils import misc

misc = misc()
from sid import download_image
# to get image provider info under variable 'dic'
from sid.config import *

#number of processors to be used
ncore = 6
list_file = '/local/AIX/tristan.harmel/project/acix/AERONETOC_Matchups_List.xlsx'
sites = pd.read_excel(list_file)  # , sep=' ')
command = []
list = []
for idx, site in sites.iterrows():
    if site.iloc[0] != site.iloc[0]:
        continue
    name, lat, lon, date_raw, time, basename = site.iloc[0:6]
    
    # get date in pratical format
    date = datetime.datetime.strptime(date_raw, '%d-%m-%Y') + datetime.timedelta(hours=time)

    sensor = misc.get_sensor(basename)
    if sensor == None:
        print('non standard image, not processed: ', basename)
        continue
    productimage = sensor[1]
    sat = sensor[2]
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
    # check if image is already downloaded
    if not os.path.exists(file):
        list.append(file)
        fromdate = date.strftime('%Y-%m-%d')
        todate = datetime.datetime.strftime(date + datetime.timedelta(days=1), '%Y-%m-%d')
        cloudmax = str(100)

        script = dic[productimage]['script']
        write = dic[productimage]['path']
        auth = dic[productimage]['auth']
        tile = misc.get_tile(basename)
        command.append([script, lat, lon, write, auth, tile, sat, cloudmax, fromdate, todate, productimage])

        # download image
        #download_image.mp_worker(command)

p = multiprocessing.Pool(ncore)
p.map(download_image.mp_worker, command)