'''
command to process images over the aeronet-oc sites
'''
import os, sys
import pandas as pd
import glob

from grs import grs_process

sitefile = sys.argv[1]
idir = '/nfs/DD/S2/L1/ESA/'
odir = sys.argv[2]  # '/nfs/DP/S2/L2/GRS/aeronet-oc/netcdf'
if not os.path.exists(odir):
    os.makedirs(odir)

sites = pd.read_csv(sitefile, sep=' ')

full_tile = True
if full_tile:
    w, h = 200, 200
else:
    w, h = 1, 1

lev = 'L2grs'
aerosol = 'cams_forecast'
fjunk = os.path.join(odir, 'list_junk_files.txt')
Nimage = 1
noclobber = True
resolution = None
aeronet_file = 'no'
aot550 = 0.1;
angstrom = 0.5

lonmin, lonmax = -180, 180
latmin, latmax = -90,-21.13
wkt_rect = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
              + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
              + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"


def wktbox(center_lon, center_lat, width=1, height=1):
    '''

    :param center_lon: decimal longitude
    :param center_lat: decimal latitude
    :param width: width of the box in km
    :param height: haight of the box in km
    :return: wkt of the box centered on provided coordinates
    '''
    from math import sqrt, atan, pi
    import pyproj
    geod = pyproj.Geod(ellps='WGS84')
    width, height = width * 1000, height * 1000
    rect_diag = sqrt(width ** 2 + height ** 2)

    azimuth1 = atan(width / height)
    azimuth2 = atan(-width / height)
    azimuth3 = atan(width / height) + pi  # first point + 180 degrees
    azimuth4 = atan(-width / height) + pi  # second point + 180 degrees

    pt1_lon, pt1_lat, _ = geod.fwd(center_lon, center_lat, azimuth1 * 180 / pi, rect_diag)
    pt2_lon, pt2_lat, _ = geod.fwd(center_lon, center_lat, azimuth2 * 180 / pi, rect_diag)
    pt3_lon, pt3_lat, _ = geod.fwd(center_lon, center_lat, azimuth3 * 180 / pi, rect_diag)
    pt4_lon, pt4_lat, _ = geod.fwd(center_lon, center_lat, azimuth4 * 180 / pi, rect_diag)

    wkt_point = 'POINT (%.6f %.6f)' % (center_lon, center_lat)
    wkt_poly = 'POLYGON (( %.6f %.6f, %.6f %.6f, %.6f %.6f, %.6f %.6f, %.6f %.6f ))' % (
        pt1_lon, pt1_lat, pt2_lon, pt2_lat, pt3_lon, pt3_lat, pt4_lon, pt4_lat, pt1_lon, pt1_lat)
    return wkt_poly


def set_ofile(file, odir='', outfile=None, level_name='l2grs', suffix=''):
    ##################################
    # File naming convention
    ##################################

    # if outfile == None:
    lev = level_name

    outfile = file.replace('L1C', lev)
    outfile = outfile.replace('.SAFE', '').rstrip('/')
    outfile = outfile.replace('.zip', '').rstrip('/')
    outfile = outfile.replace('L1TP', lev)
    outfile = outfile.replace('.txt', '').rstrip('/')

    if ('S2A' in outfile): sensor = 'S2A'
    if ('S2B' in outfile): sensor = 'S2B'
    if ('LC08' in outfile): sensor = 'LANDSAT_8'
    if ('LE07' in outfile): sensor = 'LANDSAT_7'
    if ('LT05' in outfile): sensor = 'LANDSAT_5'

    return os.path.join(odir, outfile), sensor


def chunk(it, n):
    try:
        while True:
            xs = []  # The buffer to hold the next n items
            for _ in range(n):
                xs.append(next(it))
            yield xs
    except StopIteration:
        yield xs


# load list of raw image files producing exception during grs process (causes to be investigated)
try:
    with open(fjunk) as f:
        junkfiles = f.read().splitlines()
except:
    junkfiles = []
    with open(fjunk, 'w'):
        pass

for idx, row in sites.iterrows():
    print(row)
    site = row.site
    altitude = row.alt
    imgs = glob.glob(idir + '*' + row.tile + '*')
    print(imgs.__len__())

    if False:
        wkt = wktbox(row.lon, row.lat,width=w,height=h)
    else:
        wkt = wkt_rect

    # ----------------------
    # remove file names which won't be processed whatsoever:
    imgs_tbp = []
    for file in imgs:
        # TODO remove this part to process images from 2019
        if '201901' in file:
            continue

        if file in junkfiles:
            continue

        basename = os.path.basename(file)
        if 'incomplete' in basename:
            continue

        outfile, sensor = set_ofile(basename, odir=odir)
        print(outfile, sensor)
        # if os.path.isfile(outfile + ".dim") & os.path.isdir(outfile + ".data") & noclobber:
        if os.path.isfile(outfile + ".dim") & noclobber:
            print('File ' + outfile + ' already processed; skip!')
            continue
        if os.path.isfile(outfile + ".dim.incomplete"):# & False:
            print('found incomplete File ' + outfile + '; skip!')
            continue

        imgs_tbp.append(file)
    if imgs_tbp == []:
        continue
    # ----------------------
    isuccess = 1
    for files in chunk(iter(imgs_tbp), Nimage):

        for file in files:

            # unzip image file
            unzip = False
            if os.path.splitext(file)[-1] == '.zip':
                unzip = True


            basename = os.path.basename(file)
            outfile, sensor = set_ofile(basename, odir=odir)
            print('-------------------------------')
            print('call grs for ', outfile, sensor)
            print('-------------------------------')

            #check if already partially processed, if so get startrow value
            checksum = outfile+'.checksum'
            startrow=0
            try:
                with open(checksum) as f:
                    checkdata = f.read().splitlines()
                for s in checkdata:
                    ss=s.split()
                    if ss[0]=='startrow':
                        startrow=int(ss[1])
            except:
                pass

            try:
                grs_process.process().execute(file, outfile, sensor, wkt, altitude=altitude, aerosol=aerosol,
                                              gdm=None, aeronet_file=aeronet_file, resolution=resolution,
                                              aot550=aot550, angstrom=angstrom, unzip=unzip, startrow=startrow)
                isuccess += 1
            except:
                # TODO note file name into log file
                print('-------------------------------')
                print('error for file  ', file, ' skip')
                print('-------------------------------')
                with open(fjunk, "a") as myfile:
                    myfile.write(file + '\n')
                continue

        # comment next lines if you want to process the full series of images
        # warning: this can be very (prohibitively) consuming in memory,
        # recommended use: set Nimage parameter and call this code within different subprocess
        if isuccess > Nimage:
            sys.exit()
