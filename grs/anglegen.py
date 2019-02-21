import sys
from subprocess import call
from esasnappy import ProductUtils, ProductIO

from .s2_angle import *

class angle_generator:
    def landsat(self, l2h,tartmp=None):
        '''
        Generate files for sensor and solar angles if files do not exist.
        Copy angles data into l2h.product.

        Arguments:
          * ``l2h`` -- Object for level 2 data; see "utils" module
        '''

        # set no data value
        nodatavalue = -99

        # anggen: True if angle file generated
        # (used to add files to 'image.tgz' folder
        anggen=False

        # Sun angles
        ang_type = 'solar'
        suff = ang_type + '_' + l2h.sensordata.angle_names[0] + '.img.hdr'
        ang_file = l2h.headerfile.replace('MTL.txt', suff)
        ang_header = l2h.headerfile.replace('MTL.txt', 'ANG.txt')
        if not os.path.isfile(ang_file):
            arg = ' '.join([ang_header, 'SOLAR', '1', '-b', \
                            l2h.sensordata.angle_names[0][-1],'-f',str(nodatavalue)])
            call(l2h.sensordata.angle_processor + ' ' + arg, shell=True)
            anggen=True
            # add generated file to tartmp
            if tartmp != None:
                tartmp.add(ang_file)
                tartmp.add(ang_file.replace('.hdr',''))

        #-- add to l2h.product
        ang = ProductIO.readProduct(ang_file)
        for band in ang.getBandNames():

            bandname=band.replace('Azimuth','sun_azimuth').replace('Zenith','sun_zenith')
            print("copying " + band + 'to ' + bandname)
            ProductUtils.copyBand(band, ang, bandname, l2h.product, True)
            l2h.product.getBand(bandname).setScalingFactor(0.01)
            l2h.product.getBand(bandname).setNoDataValue(nodatavalue)
            l2h.product.getBand(bandname).setNoDataValueUsed(True)

        # Band viewing angles
        ang_type = 'sensor'
        for band_name,ang_name in zip(l2h.band_names, l2h.sensordata.angle_names):
            suff = ang_type + '_' + ang_name + '.img.hdr'
            ang_file = l2h.headerfile.replace('MTL.txt', suff)
            print(ang_file)

            if not os.path.isfile(ang_file):
                arg = ' '.join([ang_header, 'SATELLITE', '1', '-b', ang_name[-1],'-f',str(nodatavalue)])
                call(l2h.sensordata.angle_processor + ' ' + arg, shell=True)
                # add generated file to tartmp

                if tartmp != None:
                    tartmp.add(ang_file)
                    tartmp.add(ang_file.replace('.hdr',''))

        #-- add to l2h.product
            ang = ProductIO.readProduct(ang_file)
            for band in ang.getBandNames():
                bandname=band+'_'+band_name
                print("copying " + band + 'to ' + bandname)
                ProductUtils.copyBand(band, ang, bandname, l2h.product, True)
                l2h.product.getBand(bandname).setScalingFactor(0.01)
                l2h.product.getBand(bandname).setNoDataValue(nodatavalue)
                l2h.product.getBand(bandname).setNoDataValueUsed(True)
        return anggen

    def landsat_tm(self, l2h):
        '''
        Angles computation for Landsat 4, 5, 7
        Generate files for sensor and solar angles if files do not exist.
        Copy angles data into l2h.product.

        Arguments:
          * ``l2h`` -- Object for level 2 data; see "utils" module

        Options:
          * ``-s 1`` -- subsampling option set to 1 (no subsampling)

        '''

        # Call angle processor executable
        nodatavalue = 0
        ang_header = l2h.headerfile.replace('MTL.txt', 'ANG.txt')
        ang_file = os.path.join(os.path.dirname(l2h.headerfile),'angle_solar_B01.img.hdr')

        for band_name,ang_name in zip(l2h.band_names, l2h.sensordata.angle_names):
            solfile = ang_file.replace('B01',ang_name)
            senfile = solfile.replace('solar','sensor')
            if not os.path.isfile(solfile) or not os.path.isfile(senfile):
                arg = ' '.join([ang_header, '-s 1'])
                call(l2h.sensordata.angle_processor + ' ' + arg, shell=True)

            #-- add to l2h.product

            # Sensor angles
            ang = ProductIO.readProduct(senfile)
            for angband in ang.getBandNames():
                bandname=angband.replace('Band_2','Azimuth').replace('Band_1','Zenith') + '_'+band_name
                print("copying " + angband + 'to ' + bandname)
                ProductUtils.copyBand(angband, ang, bandname, l2h.product, True)
                l2h.product.getBand(bandname).setScalingFactor(0.01)
                l2h.product.getBand(bandname).setNoDataValue(nodatavalue)
                l2h.product.getBand(bandname).setNoDataValueUsed(True)

        # Sun angles (keep only data for one given band; no much variations between bands
        ang = ProductIO.readProduct(solfile)
        for angband in ang.getBandNames():
            bandname = angband.replace('Band_2', 'sun_azimuth').replace('Band_1', 'sun_zenith')
            print("copying " + angband + 'to ' + bandname)
            ProductUtils.copyBand(angband, ang, bandname, l2h.product, True)
            l2h.product.getBand(bandname).setScalingFactor(0.01)
            l2h.product.getBand(bandname).setNoDataValue(nodatavalue)
            l2h.product.getBand(bandname).setNoDataValueUsed(True)

    def sentinel2(self, l2h):
        ''' Generate files for sensor angles if files do not exist.
        Copy angles data into l2h.product.


        Arguments:
          * ``l2h`` -- Object for level 2 data; see "utils" module '''


        import glob

        ##################################
        # GENERATE BAND ANGLES
        ##################################
        # get xml file of tile to be processed
        root = os.path.abspath(str(l2h.product.getFileLocation()))
        root = os.path.join(root, 'GRANULE')
        granule = glob.glob1(root, '*L1C*')
        if granule.__len__() > 1:
            print('STOP! should process single tile but several tiles found')
            sys.exit()
        root = os.path.join(root, granule[0])
        XML_File = os.path.join(root, glob.glob1(root, '*MTD*xml')[0])

        outdir = os.path.join(root, 'ANG_DATA')
        print(outdir)
        if not os.path.exists(outdir) or os.listdir(outdir) == "":
            print('creating angle files for ', outdir)
            s2angle().angle_writer(XML_File)

        ##################################
        # ADD BAND ANGLES TO PRODUCT
        ##################################
        angroot = os.path.join(root, 'ANG_DATA')
        angfiles = glob.glob1(angroot, '*hdr')

        for angfile in angfiles:
            bandname = angfile[-7:-4].replace('B0', 'B')
            ang = ProductIO.readProduct(os.path.join(angroot, angfile))
            # ProductUtils.copyGeoCoding(product, ang)
            for band in ang.getBandNames():
                print("copying " + band + '_' + bandname)
                ProductUtils.copyBand(band, ang, band + '_' + bandname, l2h.product, True)


    def add_band(self,):
        print()
