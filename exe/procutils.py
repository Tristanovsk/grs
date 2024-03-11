'''utils module dedicated to processing of massive dataset'''

import os, sys
import re
import numpy as np
import pandas as pd
import glob
import datetime


class misc:
    '''
    Miscellaneous utilities
    '''
    def __init__(self):
        pass

    def wktbox(self, center_lon, center_lat, width=1, height=1):
        '''

        :param center_lon: decimal longitude
        :param center_lat: decimal latitude
        :param width: width of the box in km
        :param height: height of the box in km
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

    def get_sensor(self, file):
        '''
        Get sensor type from file name
        :param file: file in standard naming
        :return: sensor type
        '''
        file = os.path.basename(file)

        if ('S2A' in file):
            sensor = ('S2A', 'S2_ESA', None)
        elif ('S2B' in file):
            sensor = ('S2B', 'S2_ESA', None)
        elif ('LC08' in file) | bool(re.search(r"L[C,O]8",file)):
            sensor = ('LANDSAT_8', 'Landsat_USGS', 'LC8')
        elif ('LE07' in file) | ('LE7' in file):
            sensor = ('LANDSAT_7', 'Landsat_USGS', 'LE7')
        elif ('LT05' in file) | ('LT5' in file):
            sensor = ('LANDSAT_5', 'Landsat_USGS', 'LT5')
        # TODO add to log file
        else:
            print('sensor not recognized from input file')
            sensor = None
            # sys.exit(-1)
        return sensor

    def get_tile(self, file):
        '''
        Get tile from file name
        :param file: file in standard naming
        :return: tile ID
        '''
        file = os.path.basename(file)
        sensor = ''
        if ('S2A' in file) | ('S2B' in file):
            tile = file.split('_')[5][-5:]
        elif bool(re.search(r"L[C,O]8",file)) | ('LE7' in file) | ('LT5' in file):
            tile = file[3:9]
        elif ('LC08' in file) | ('LE07' in file) | ('LT05' in file):
            tile = file[4:10]
        # TODO add to log file
        else:
            print('sensor not recognized from input file')
            # sys.exit(-1)
        return tile

    def set_ofile(self, file, odir='', level_name='L2GRS', suffix=''):
        ''' get satellite type andset output file name'''
        ##################################
        # File naming convention
        ##################################

        lev = level_name

        outfile = file.replace('L1C', lev)
        outfile = outfile.replace('L1GT', lev)
        outfile = outfile.replace('L1TP', lev)
        # remove extension
        outfile = os.path.splitext(outfile)[0]
        # outfile = outfile.replace('.SAFE', '').rstrip('/')
        # outfile = outfile.replace('.zip', '').rstrip('/')
        # outfile = outfile.replace('.txt', '').rstrip('/')
        #self.get_tile(file),
        path = os.path.join(odir, outfile + suffix)
        print(path)
        return path

    def chunk(self, it, n):
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

class multi_process:
    '''
    Utilities for multicore processing
    '''
    def __init__(self):
        pass

    def grs_call(self,p):
        args,fjunk = p
        for arg in args:
            file_tbp, outfile, wkt, altitude, aerosol, aeronet_file, ancillary, resolution, \
            aot550, angstrom, mem_safe, unzip, untar, startrow, allpixels, angleonly = arg
            print('yop',file_tbp)
            #return

            try:
                from grs import grs_process
                grs_process.Process().execute(file_tbp, outfile, wkt, altitude=altitude, aerosol=aerosol, ancillary=ancillary,
                                              dem=True, aeronet_file=aeronet_file, resolution=resolution,
                                              aot550=aot550, angstrom=angstrom, memory_safe=mem_safe, unzip=unzip, untar=untar,
                                              startrow=startrow, allpixels=allpixels, angleonly=angleonly)
            except:
                print('-------------------------------')
                print('error for file  ', file_tbp, ' skip')
                print('-------------------------------')
                with open(fjunk, "a") as myfile:
                    myfile.write(file_tbp + ' error during grs \n')
                continue
        # here sys.exit instead of "return" to terminate and close snappy and free memory
        sys.exit()
        return

    def grs_cnes(self,arg):

        # for arg in args:
        print('arg', arg)
        file_tbp, outfile, aerosol, aeronet_file, ancillary, resolution, \
        dem, maja_xml, waterdetect_file, \
        aot550, angstrom, mem_safe, allpixels, angleonly = arg
        print('start process of ',file_tbp)

        from grs import grs_process

        grs_process.Process().execute(file_tbp, outfile, aerosol=aerosol, ancillary=ancillary,
                                      dem=dem, aeronet_file=aeronet_file, resolution=resolution,
                                      maja_xml=maja_xml, waterdetect_file=waterdetect_file,
                                      aot550=aot550, angstrom=angstrom, memory_safe=mem_safe,
                                      allpixels=allpixels, angleonly=angleonly)

        # here sys.exit instead of "return" to terminate and close snappy and free memory
        sys.exit()
        return

