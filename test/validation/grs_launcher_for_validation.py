'''
To provide file to run a specific grs version on validation dataset:
```
python test/validation/grs_launcher_for_validation.py > run_test.sh
```
'''

import os
import glob
import xarray as xr
import pandas as pd

from datetime import datetime, timedelta
import datetime as dt
import grs

#print(f'-grs: {grs.__version__}')
opj = os.path.join


l1c_subset_dir='/work/datalake/watcal/validation/L1C/'
l2a_dir='/work/datalake/watcal/validation/L2A/'

cams_dir = '/work/datalake/static_aux/CAMS/OBS2CO_GRS_v2'

dem_dir = '/work/datalake/static_aux/MNT/COP-DEM_GLO-30-DGED_S2_tiles'
start_date='2019-01-01'
end_date='2020-12-31'
site= 'BEFR'
resolution=20
cc_max = 0.6
no_clobber = True

# TODO add DEM file
#dem_file = opj(dem_dir,'COP-DEM_GLO-30-DGED_'+tile+'.tif')

# --------------------
# set options
options = ''
if no_clobber:
    options = options + ' --no_clobber'
print('#!/bin/bash')
for year_ in range(2020, 2023):
    year = str(year_)

    files = glob.glob(opj(l1c_subset_dir, site, year, '*.nc'))
    files.sort()
    files
    for file in files:

        basename = os.path.basename(file).split('.')[0].replace('L1C', 'L2A')
        info = basename.split('_')
        time = dt.datetime.strptime(info[2], '%Y%m%dT%H%M%S')
        date_str = time.strftime('%Y-%m-%d')
        date_rep = time.strftime('%Y/%m/%d')
        tile = info[5][1:]
        l2a_dir_ = opj(l2a_dir, site, date_rep)
        l2a_filename = basename + '_V' + grs.__version__
        l2a_path = opj(l2a_dir_, l2a_filename)

        if os.path.exists(l2a_path):
            #    pass
            continue


        cams_file = opj(cams_dir, date_rep, date_str+'-cams-global-atmospheric-composition-forecasts.nc')
        # cams_file = '/home/harmel/Downloads/cams_forecast_2024-09-08_beta.nc'


        print(f'grs {file} -o {l2a_dir_}' +  # f'--dem_file {dem_file} '+
              f' --cams_file {cams_file}  --max_cloud_cover {cc_max} --resolution {resolution}' + options)