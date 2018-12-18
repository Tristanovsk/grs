'''
command to process images over the aeronet-oc sites
'''

import grs
import pandas as pd
import glob

sitefile='exe/aeronet-oc_sites.txt'
sites=pd.read_csv(sitefile)

idir='/nfs/HYAX/imagerie/S2/L1/ESA/'

for idx, row in sites.iterrows():
    print(row.site,row.lon)
    glob.glob(idir+'*'+row.tile)

for x in sites:
    print(x)
glob.glob(idir+'*'+tile)



sensor='S2A'
odir='/net/axsimagerie/mnt/datas/imagerie/S2/L2/GRS/'
lev='L2grs_test'
aerosol='cams_forecast'

noclobber = False
resolution = None
aeronet_file = 'no'
aot550=0.1; angstrom=0.5

altitude=0
lonmax = 180
lonmin = -180
latmax = 90
latmin = -90



grs.grs_process.process().execute(file, outfile, sensor, wkt, altitude=altitude, aerosol=aerosol,
                      gdm=None, aeronet_file=aeronet_file, resolution=resolution,
                      aot550=args['--aot550'], angstrom=args['--angstrom'])
