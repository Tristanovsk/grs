'''
command to load L1C and L2A images on datalake from the given list

example:
python3 call_amalthee.py list_grs_cnes_template.csv
'''

import os, sys
import numpy as np
import pandas as pd
import glob
from datetime import datetime, timedelta
import time

# CNES lib for datalake managment
# import lxml

from libamalthee import Amalthee

L2A = Amalthee('theia')
L1C = Amalthee('peps')

# start_date, end_date = '2018-01-01', '2018-01-20'
# tile, lon, lat = '33PVR', '14.6', '14'

sitefile = 'exe/list_grs_cnes_obs2mod.csv'
sitefile = sys.argv[1]
sites = pd.read_csv(sitefile)

# --------------------------------------------------------------------------------

for idx, site in sites.iterrows():
    # load row of list file
    if site.iloc[0] == 0:
        continue
    name, start_date, end_date, sat, tile, lat, lon, resolution = site.iloc[1:]
    if start_date == end_date:
        end_date = (datetime.strptime(end_date, '%Y-%m-%d').date() + timedelta(days=1)).__str__()

    parameters = {"productType": "S2MSI1C", "tileid": tile}
    L1C.search("S2ST", start_date, end_date, parameters)
    L1C.fill_datalake()

    parameters = {'processingLevel': 'LEVEL2A', 'lon': str(lon), 'lat': str(lat)}
    L2A.search("SENTINEL2", start_date, end_date, parameters)
    L2A.fill_datalake()

    finished = False
    while not finished:
        L1C.check_datalake()
        L2A.check_datalake()
        L2A.products.loc[L2A.products.state == 'failed','available']=True
        finished = all(L2A.products.available) and all(L1C.products.available)
        time.sleep(3)
    print('Finished!!')
