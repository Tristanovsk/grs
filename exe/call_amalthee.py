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


# start_date, end_date = '2021-03-01', '2021-03-30'
# tile, lon, lat = '31TGK', '14.6', '14'

sitefile = 'exe/list/list_grs_cnes_seine.csv' # 'exe/list_grs_cnes_seine.csv' #
sitefile = sys.argv[1]
sites = pd.read_csv(sitefile)

# --------------------------------------------------------------------------------

for idx, site in sites.iterrows():
    # load row of list file
    if site.iloc[0] == 0:
        continue

    L2A = Amalthee('theia')
    L1C = Amalthee('peps')
    name, start_date, end_date, sat, tile, resolution, flag = site.iloc[1:8]
    if start_date == end_date:
        end_date = (datetime.strptime(end_date, '%Y-%m-%d').date() + timedelta(days=1)).__str__()

    parameters = {"productType": "S2MSI1C", "tileid": tile}
    L1C.search("S2ST", start_date, end_date, parameters)
    idL1C = L1C.fill_datalake()
    L1C.check_datalake()
    parameters = {'processingLevel': 'LEVEL2A', 'location':'T'+tile} #'lon': str(lon), 'lat': str(lat)}
    L2A.search("SENTINEL2", start_date, end_date, parameters)
    idL2A = L2A.fill_datalake()

    # if idx == 1:
    #     break
    finished = False
    finished_L1C = False
    finished_L2A = False
    iwait = 0
    while not finished:
        print('tile',tile)
        res = L1C.check_datalake()
        print('L1C',res)
        try:
            if (res['status'] == 'no_request_made') or (res['status'] == 'done'):
                finished_L1C = True
                print('L1C job finished or canceled', res)
                L1C.delete_request(idL1C)
        except:
            pass
        res = L2A.check_datalake()
        print('L2A',res)

        try:
            if (res['status'] == 'no_request_made') or (res['status'] == 'done'):
                finished_L2A = True
                print('L2A job finished or canceled', res)
                L2A.delete_request(idL2A)
        except:
            pass
        # L2A.products.loc[L2A.products.state == 'failed', 'available'] = True
        #finished = (all(L2A.products.available) or finished_L2A) and (all(L1C.products.available) or finished_L1C)
        finished = ( finished_L2A) and ( finished_L1C)

        if finished:
            break

        time.sleep(123)

        iwait += 1
        if iwait > 20:
            print('time limit exceeded')
            finished = True

    print('Finished!!')
