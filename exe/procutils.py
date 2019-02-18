'''utils module dedicated to processing of massive dataset'''

import os, sys
import numpy as np
import pandas as pd
import glob
import datetime


class misc:
    def __init__(self):
        pass

    def wktbox(self,center_lon, center_lat, width=1, height=1):
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


    def get_sensor(self,file):
        '''
        Get sensor type from file name
        :param file: file in standard naming
        :return: sensor type
        '''
        file = os.path.basename(file)
        sensor=''
        if ('S2A' in file): sensor = 'S2A'
        elif ('S2B' in file): sensor = 'S2B'
        elif ('LC08' in file) | ('LC8' in file): sensor = 'LANDSAT_8'
        elif ('LE07' in file): sensor = 'LANDSAT_7'
        elif ('LT05' in file): sensor = 'LANDSAT_5'
        # TODO add to log file
        else:
            print('sensor not recognized from input file')
            #sys.exit(-1)
        return sensor

    def set_ofile(self,file, odir='', outfile=None, level_name='l2grs', suffix=''):
        ''' get satellite type andset output file name'''
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
        print(outfile)
        sensor=''
        if ('S2A' in outfile): sensor = 'S2A'
        if ('S2B' in outfile): sensor = 'S2B'
        if ('LC08' in outfile) | ('LC8' in outfile) : sensor = 'LANDSAT_8'
        if ('LE07' in outfile): sensor = 'LANDSAT_7'
        if ('LT05' in outfile): sensor = 'LANDSAT_5'

        return os.path.join(odir, outfile), sensor


    def chunk(self,it, n):
        '''
        return a tuple of n items found in it
        :param it:
        :param n:
        :return:
        '''
        try:
            while True:
                xs = []  # The buffer to hold the next n items
                for _ in range(n):
                    xs.append(next(it))
                yield xs
        except StopIteration:
            yield xs