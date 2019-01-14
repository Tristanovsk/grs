'''
command to process images over the aeronet-oc sites
'''
import os,sys
import pandas as pd
import glob
from grs import grs_process

sitefile='exe/aeronet-oc_sites.txt'
sites=pd.read_csv(sitefile,sep=' ')

idir='/nfs/DD/S2/L1/ESA/'
odir='/nfs/DP/S2/L2/GRS/aeronet-oc'
lev='L2grs'
aerosol='cams_forecast'

noclobber = True
resolution = None
aeronet_file = 'no'
aot550=0.1; angstrom=0.5

def wktbox(center_lon, center_lat, width=500,height=500):
    '''

    :param center_lon: decimal longitude
    :param center_lat: decimal latitude
    :param width: width of the box in m
    :param height: haight of the box in m
    :return: wkt of the box centered on provided coordinates
    '''
    from math import sqrt,atan,pi
    import pyproj
    geod = pyproj.Geod(ellps='WGS84')

    rect_diag = sqrt( width**2 + height**2 )

    azimuth1 = atan(width/height)
    azimuth2 = atan(-width/height)
    azimuth3 = atan(width/height)+pi # first point + 180 degrees
    azimuth4 = atan(-width/height)+pi # second point + 180 degrees

    pt1_lon, pt1_lat, _ = geod.fwd(center_lon, center_lat, azimuth1*180/pi, rect_diag)
    pt2_lon, pt2_lat, _ = geod.fwd(center_lon, center_lat, azimuth2*180/pi, rect_diag)
    pt3_lon, pt3_lat, _ = geod.fwd(center_lon, center_lat, azimuth3*180/pi, rect_diag)
    pt4_lon, pt4_lat, _ = geod.fwd(center_lon, center_lat, azimuth4*180/pi, rect_diag)

    wkt_point = 'POINT (%.6f %.6f)' % (center_lon, center_lat)
    wkt_poly = 'POLYGON (( %.6f %.6f, %.6f %.6f, %.6f %.6f, %.6f %.6f, %.6f %.6f ))' % (pt1_lon, pt1_lat, pt2_lon, pt2_lat, pt3_lon, pt3_lat, pt4_lon, pt4_lat, pt1_lon, pt1_lat)
    return wkt_poly

def set_ofile(file,odir='',outfile=None,level_name='l2grs',suffix=''):
    ##################################
    # File naming convention
    ##################################

    #if outfile == None:
    lev = level_name

    outfile = file.replace('L1C', lev)
    outfile = outfile.replace('.SAFE', '').rstrip('/')
    outfile = outfile.replace('.zip', '').rstrip('/')
    outfile = outfile.replace('L1TP', lev)
    outfile = outfile.replace('.txt', '').rstrip('/')

    if('S2A'in outfile):sensor='S2A'
    if('S2B'in outfile):sensor='S2B'
    if('LC08'in outfile):sensor='LANDSAT_8'
    if('LE07'in outfile):sensor='LANDSAT_7'
    if('LT05'in outfile):sensor='LANDSAT_5'

    return os.path.join(odir, outfile),sensor

for idx, row in sites.iterrows():
    site=row.site
    altitude=row.alt
    imgs=glob.glob(idir+'*'+row.tile+'*')
    print(imgs.__len__())
    wkt=wktbox(row.lon,row.lat)

    for file in imgs:
# TODO remove this part to process images from 2019
        if '201901' in file:
            continue

        basename = os.path.basename(file)
        if 'incomplete' in basename:
            continue
        outfile,sensor=set_ofile(basename,odir=odir)
        print(outfile,sensor)
        if os.path.isfile(outfile + ".dim") & os.path.isdir(outfile + ".data") & noclobber:
                print('File ' + outfile + ' already processed; skip!')
                continue

        #unzip image file
        unzip=False
        if os.path.splitext(file)[-1] == '.zip':
            unzip=True

        grs_process.process().execute(file, outfile, sensor, wkt, altitude=altitude, aerosol=aerosol,
                                          gdm=None, aeronet_file=aeronet_file, resolution=resolution,
                                          aot550=aot550, angstrom=angstrom, unzip=unzip)
