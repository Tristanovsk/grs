'''
Main program
'''

from pathlib import Path
from esasnappy import ProductData, ProductIO
import logging
from logging.handlers import RotatingFileHandler

import os, shutil
import zipfile
import tarfile
import glob
import numpy as np
import xarray as xr

from . import config as cfg
from . import acutils
from . import auxdata
from . import utils
from .anglegen import *
from .fortran.grs import main_algo as grs_solver
from .fortran.grs_a import main_algo as grs_a_solver


class process:
    ''' '''

    def __init__(self):
        pass

    def execute(self, file, outfile, wkt=None, sensor=None, aerosol='default', ancillary=None, altitude=0,
                dem=True, aeronet_file=None, aot550=0.1, angstrom=1, resolution=None, unzip=False, untar=False,
                startrow=0, allpixels=False, maja_xml=None, waterdetect_file=None, waterdetect_only=False,
                memory_safe=False, angleonly=False, grs_a=False, output='Rrs', logfile="log.txt", log_level="INFO", 
                xblock=512, yblock=512):
        '''
        Main program calling all GRS steps

        :param file: Input file to be processed
        :param outfile: Absolute path of the output file
        :param wkt: Well-Known-Text format defining the area of interest for which the image is subset
        :param sensor: Set the sensor type: S2A, S2B, LANDSAT_5, LANDSAT_7, LANDSAT_8
                    (by default sensor type is retrieved from input file name)
        :param aerosol: aerosol data base to use within the processing
                   DB: cams_forecast, cams_reanalysis, cds_forecast, aeronet, user_model, default
        :param ancillary: if None, value is set to that of aerosol
        :param altitude: provide altitude if `dem` is set as `False`
        :param dem: if True digital elevation model is applied for per-pixel pressure calculation (data from SNAP/SRTM)
        :param aeronet_file: optional aeronet file to be used for aerosol calculations
        :param maja_xml: optional use of mask from MAJA L2A images, path to xml ID of the L2A image
        :param waterdetect_file: optional use of water mask from waterdetect algorithm,
                    path to the appropriate WaterDetect data file
        :param waterdetect_only: if True and waterdetect file is provided, process only the pixels masked as "water"
        :param resolution: pixel resolution in meters (integer)
        :param unzip: if True input file is unzipped before processing,
                      NB: unzipped files are removed at the end of the process
        :param startrow: row number of the resampled and subset image on which the process starts, recommended value 0
                        NB: this option is used to in the context of operational processing of massive dataset
        :param allpixels: force to process all pixels even they are flagged as "Vegetation" or "Non-water"
        :param angleonly: if true, grs is used to compute angle parameters only (no atmo correction is applied)
        :param output: set the unit of the retrievals:

                 * 'Lwn', normalized water-leaving radiance (in  :math:`mW cm^{-2} sr^{-1} \mu m^{-1})`

                 * 'Rrs', remote sensing reflectance (in  :math:`sr^{-1}`)

                 {default: 'Rrs'}
        :param grs_a: switch to grs-a algorithm (Lwn and aerosol) if True

        :return:
        '''


        ##################################
        # Get sensor auxiliary data
        ##################################

        logging.info('Get sensor auxiliary data')
        _utils = utils.utils()
        if sensor == None:
            sensor = _utils.get_sensor(file)
        sensordata = auxdata.sensordata(sensor)
        if resolution == None:
            resolution = sensordata.resolution
        indband = sensordata.indband

        if ancillary == None:
            ancillary = aerosol

        ##################################
        # Read L1C product
        ##################################


        file_orig = file
        # unzip if needed
        if unzip:
            logging.info('unzipping...')
            tmpzip = zipfile.ZipFile(file)
            tmpzip.extractall(cfg.tmp_dir)
            file = os.path.join(cfg.tmp_dir, tmpzip.namelist()[0])
        tartmp = None
        if untar:
            basename = os.path.basename(file).replace('.tgz', '')
            basename = os.path.basename(basename).replace('.tar.gz', '')
            tmp_dir = os.path.join(cfg.tmp_dir, basename)
            # open tar archive to extract files for data loading
            tmpzip = tarfile.open(file)
            tmpzip.extractall(tmp_dir)
            file = glob.glob(os.path.join(tmp_dir, '*MTL.*'))[0]
            # open tar archive to add potential file (e.g., angle files) - InvalidHeaderError
            if not any(['solar' in f for f in glob.glob(os.path.join(tmp_dir, '*'))]):
                tartmp = tarfile.open(os.path.join(cfg.tmp_dir, os.path.basename(file_orig)), 'w:gz')

        logging.info("Reading...")
        logging.info(file)
        product = ProductIO.readProduct(file)

        ##################################
        # Generate l2h object
        ##################################
        l2h = utils.info(product, sensordata, aerosol, ancillary, output)
        l2h.headerfile = file

        ##################################
        # GET METADATA
        ##################################
        logging.info('getting metadata...')
        # TODO clean up this part and other metadata to be loaded
        if 'S2' in sensor:
            meta = l2h.product.getMetadataRoot().getElement('Level-1C_User_Product').getElement(
                'General_Info').getElement(
                'Product_Image_Characteristics').getElement('Reflectance_Conversion')
            l2h.U = float(str(meta.getAttribute('U').getData()))
            l2h.solar_irr = np.zeros(len(indband), dtype=np.float32)
            for i, iband in zip(range(len(indband)), indband):
                l2h.solar_irr[i] = float(str(meta.getElement('Solar_Irradiance_List').getAttributeAt(iband).getData()))

        else:
            meta = l2h.product.getMetadataRoot().getElement("L1_METADATA_FILE").getElement("IMAGE_ATTRIBUTES")
            l2h.U = float(str(meta.getAttribute('EARTH_SUN_DISTANCE').getData())) ** 2
            l2h.solar_irr = np.array(l2h.sensordata.solar_irr)[indband]

        # convert into mW cm-2 um-1
        l2h.solar_irr = l2h.solar_irr / 10

        ##################################
        # GENERATE BAND ANGLES (LANDSAT)
        ##################################
        anggen = False
        if 'LANDSAT_8' in sensor:
            anggen = angle_generator().landsat(l2h)
        elif 'LANDSAT' in sensor:
            anggen = angle_generator().landsat_tm(l2h)
        #logging.info('anggen = {}'.format(anggen))
        if anggen:
            if tartmp:
                logging.info('writing angles to input file: ' + file_orig)
                # TODO finalize this part to add angle files to original tar.gz image (e.g., LC8*.tgz)
                # copy tgz image and add angle files
                shutil.move(file_orig, os.path.join(os.path.dirname(file_orig), 'saves', os.path.basename(file_orig)))
                for filename in os.listdir(tmp_dir):
                    tartmp.add(os.path.join(tmp_dir, filename), filename)
                tartmp.close()
                shutil.move(os.path.join(cfg.tmp_dir, os.path.basename(file_orig)), file_orig)

        # stop process for landsat angle computation only
        if angleonly:
            return

        ##################################
        # RESAMPLE TO A UNIQUE RESOLUTION
        ##################################
        logging.info('resampling...')
        if 'S2' in sensor:
            if memory_safe:
                l2h.product = _utils.generic_resampler(l2h.product, resolution=resolution)  # , method='Nearest')
            else:
                l2h.product = _utils.s2_resampler(l2h.product, resolution=resolution)
        else:
            l2h.product = _utils.resampler(l2h.product, resolution=resolution)  # , upmethod='Nearest')


        ##################################
        # SUBSET TO AREA OF INTEREST
        ##################################
        logging.info('subsetting...')
        try:
            if wkt is not None:
                l2h.product = _utils.get_subset(l2h.product, wkt)
        except:
            if unzip:
                # remove unzipped files (Sentinel files)
                shutil.rmtree(file, ignore_errors=True)
            if untar:
                # remove untared files (Landsat files)
                shutil.rmtree(tmp_dir, ignore_errors=True)
            raise NameError('No data available for requested area')

        l2h.get_product_info()
        l2h.set_outfile(outfile)
        l2h.wkt, lonmin, lonmax, latmin, latmax = _utils.get_extent(l2h.product)
        l2h.crs = str(l2h.product.getBand(l2h.band_names[0]).getGeoCoding().getImageCRS())

        ##################################
        # Fetch optional mask products
        # resample for common resolution
        # subset to ROI
        ##################################
        logging.info('fetching flags...')
        maja, waterdetect = None, None
        if maja_xml:
            try:
                maja = ProductIO.readProduct(maja_xml)
                maja = _utils.resampler(maja, resolution=resolution)
                maja = _utils.get_subset(maja, wkt)

            except:
                logging.info('!!! issues with ' + maja_xml + '; please check if file exists')
                raise

        if waterdetect_file:
            try:
                waterdetect = ProductIO.readProduct(waterdetect_file)
                waterdetect = _utils.get_subset(waterdetect, wkt)
                waterdetect = _utils.resampler(waterdetect, resolution=resolution)

            except:
                logging.info('!!! issues with ' + waterdetect_file + '; please check if file exists')
                raise

        ##################################
        ## ADD ELEVATION BAND
        ##################################
        logging.info('adding elevation band...')
        if dem:
            logging.info('add elevation band')
            high_latitude = (latmax >= 60) | (latmin <= -60)
            l2h.get_elevation(high_latitude)

        else:
            l2h.elevation = np.zeros([l2h.height, l2h.width])

        ##################################
        # GET IMAGE AND RASTER PROPERTIES
        ##################################
        logging.info('load raster data...')
        l2h.get_bands(l2h.band_names)
        l2h.print_info()

        ##################################
        # SET NEW BAND FOR MASKING
        ##################################
        logging.info('add ndwi mask')
        l2h.product.addBand('ndwi', ProductData.TYPE_FLOAT32)
        l2h.ndwi_band = l2h.product.getBand('ndwi')
        l2h.ndwi_band.ensureRasterData()
        l2h.ndwi_band.loadRasterData()

        ##################################
        # GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
        ##################################
        logging.info('getting CAMS data...')
        l2h.aux = auxdata.cams()

        if ancillary != 'default':
            if l2h.aerosol == 'cds_forecast':
                cams_file=os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'),l2h.date.strftime('%m'),l2h.date.strftime('%d'),
                                 l2h.date.strftime('%Y-%m-%d') + '-cams-global-atmospheric-composition-forecasts.nc')
                if(not os.path.exists(cams_file)):
                    cams_file = os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'), l2h.date.strftime('%Y-%m') +
                                      '_month_cams-global-atmospheric-composition-forecasts.nc')
                logging.info(str(cams_file)+" cams file would be used")
                l2h.aux.get_cams_ancillary(cams_file, l2h.date, l2h.wkt, param=['msl', 'gtco3', 'tcwv', 'tcno2', 't2m'])
            else:
                cams_file = Path(os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'),
                                           l2h.date.strftime('%Y-%m') + '_month_' + l2h.ancillary + '.nc'))
                logging.info(str(cams_file)+" cams file would be used")
                # do not load here since already implemented elsewhere in CNES HPC
                # l2h.aux.load_cams_data(target, l2h.date, data_type=l2h.ancillary)
                l2h.aux.get_cams_ancillary(cams_file, l2h.date, l2h.wkt)

        ## uncomment this part to use ecmwf files provided in the .SAFE format
        # if 'S2' in sensor:
        #     l2h.aux.get_tile_dir(file)
        #     l2h.aux.get_aux_dir()
        #     l2h.aux.get_ecmwf_data()

        # get pressure at the scene altitude
        l2h.pressure_msl = l2h.aux.msl  # acutils.misc.get_pressure(altitude, l2h.aux.msl)
        if dem:
            altitude = l2h.elevation
            altitude[altitude < -200] = 0
        l2h.pressure = acutils.misc.get_pressure(altitude, l2h.pressure_msl)  # l2h.aux.pressure = l2h.pressure

        ######################################
        #      Create output l2 product
        #          'l2_product'
        ######################################
        logging.info('creating L2 output product')
        l2h.create_product(maja=maja, waterdetect=waterdetect)
        try:
            l2h.load_data()
        except:
            if unzip:
                # remove unzipped files (Sentinel files)
                shutil.rmtree(file, ignore_errors=True)
            if untar:
                # remove untared files (Landsat files)
                shutil.rmtree(tmp_dir, ignore_errors=True)
            raise NameError('No data available for requested area')
        l2h.load_flags()

        #---------
        # delete Cache
        # TODO need to find a way to delete the processed file only (keep the others)
        #l2h.deleteCache()

        #####################################
        # LOAD LUT FOR ATMOSPHERIC CORRECTION
        #####################################
        logging.info('loading lut...'+ l2h.lutfine)
        lutf = acutils.lut(l2h.band_names)
        lutc = acutils.lut(l2h.band_names)
        lutf.load_lut(l2h.lutfine, indband)
        lutc.load_lut(l2h.lutcoarse, indband)

        # reproject lut array on the angles of the image
        # angles are rounded to reduce the dims of interpolated LUT
        sza_ = _utils.remove_na(np.unique(l2h.sza.round(1)))
        vza_ = _utils.remove_na(np.unique(l2h.vza.round(1)))
        azi_ = _utils.remove_na(np.unique(l2h.razi.round(0)))
        lutf.interp_n_slice(sza_,vza_,azi_)
        lutc.interp_n_slice(sza_,vza_,azi_)
        aotlut = np.array(lutf.aot, dtype=l2h.type)


        ##################################
        # GET ANCILLARY DATA (AEROSOL)
        ##################################
        aero = acutils.aerosol()
        aot550rast = np.zeros([l2h.height, l2h.width], dtype=l2h.type, order='F')
        aotscarast = np.zeros([l2h.height, l2h.width], dtype=l2h.type, order='F')
        # ssarast = np.zeros([l2h.N,l2h.width, l2h.height], dtype=l2h.type)
        aotrast = np.zeros([l2h.N, l2h.height, l2h.width], dtype=l2h.type, order='F')
        fcoefrast = np.zeros([l2h.height, l2h.width], dtype=l2h.type, order='F')
        # set mean SSA value to adjust scattering AOT when info is not available
        ssacoef = 0.99

        if l2h.aerosol == 'cds_forecast':

            l2h.aux.get_xr_cams_cds_aerosol(cams_file, l2h, lutf, lutc)
            aot550rast = l2h.aux.aot_sca_550  # .T
            fcoefrast = l2h.aux.fcoef
            aotrast = l2h.aux.aot_grs
            aotscarast = l2h.aux.aot_sca_grs
            logging.info(f'aot550rast shape {aot550rast.shape}')

        else:
            # AERONET data
            if (l2h.aerosol == 'aeronet'):
                l2h.set_aeronetfile(aeronet_file)
                try:
                    l2h.aux.Aeronet.import_aeronet_data(aero, l2h.aeronetfile, l2h.date)
                except:
                    logging.info('Error: No aeronet data in the +/- 2day time window.')
                    sys.exit()

                l2h.aux.aot_wl = aero.wavelengths
                l2h.aux.aot = aero.aot
                aotscarast = ssacoef*l2h.aux.aot
                l2h.aux.aot550 = aero.aot550
                aot550rast.fill(l2h.aux.aot550)

            # CAMS dataset
            elif (l2h.aerosol == 'cams_forecast') | (l2h.aerosol == 'cams_reanalysis'):

                l2h.aux.get_xr_cams_aerosol(cams_file, l2h.product)
                aotscarast = ssacoef*l2h.aux.aot
                aot550rast = l2h.aux.aot550rast  # .T
                logging.info(f'aot550rast shape {aot550rast.shape}')
            # CAMS new cds dataset (available from 26 June 2018 12UTC)

            elif (l2h.aerosol == 'user_model'):
                l2h.aux.aot550 = aot550
                l2h.angstrom = angstrom
                l2h.aux.aot = l2h.aux.aot550 * (np.array(l2h.wl) / 550) ** (-l2h.angstrom)
                aotscarast = ssacoef*l2h.aux.aot
                l2h.aux.aot_wl = l2h.wl
                aot550rast.fill(l2h.aux.aot550)

            else:
                l2h.aux.aot550 = 0.1
                l2h.angstrom = 1
                l2h.aux.aot = l2h.aux.aot550 * (np.array(l2h.wl) / 550) ** (-l2h.angstrom)
                aotscarast = ssacoef*l2h.aux.aot
                l2h.aux.aot_wl = l2h.wl
                aot550rast.fill(l2h.aux.aot550)

                logging.info("No aerosol data provided, set to default: aot550=01, angstrom=1")

            # set spectral aot for satellite bands
            aero.fit_spectral_aot(l2h.aux.aot_wl, l2h.aux.aot)
            l2h.aot = aero.get_spectral_aot(np.array(l2h.wl))
            l2h.aot550 = l2h.aux.aot550

            # normalization of Cext to get spectral dependence of fine and coarse modes
            nCext_f = lutf.Cext / lutf.Cext550
            nCext_c = lutc.Cext / lutc.Cext550
            logging.info('param aerosol {nCext_f}, {nCext_c}, {l2h.aot}')
            aero.fit_aero(nCext_f, nCext_c, l2h.aot / l2h.aot550)
            logging.info(f'{aero.fcoef} {aero.fcoef.astype(l2h.type)}')
            l2h.fcoef = aero.fcoef.astype(l2h.type)
            fcoefrast.fill(l2h.fcoef[0])
            for i in range(l2h.N):
                aotrast[i].fill(l2h.aot[i])

        l2h.rot = l2h.sensordata.rot

        ####################################
        #     Set SMAC parameters for
        #    absorbing gases correction
        ####################################
        logging.info('loading SMAC algorithm...')
        smac = acutils.smac(l2h.sensordata.smac_bands, l2h.sensordata.smac_dir)
        smac.set_gas_param()
        smac.set_values(o3du=l2h.aux.o3du, h2o=l2h.aux.h2o)
        smac.set_standard_values(l2h.pressure_msl)
        l2h.aux.no2 = smac.uno2


        ######################################
        # arrays allocation
        # reshaping for fortran binding
        ######################################

        aot550guess = np.zeros(l2h.width, dtype=l2h.type)
        rtoaf = np.zeros((lutf.aot.__len__(), l2h.N, l2h.width), dtype=l2h.type, order='F')
        rtoac = np.zeros((lutc.aot.__len__(), l2h.N, l2h.width), dtype=l2h.type, order='F')
        maskpixels_ = np.full(l2h.width, 1, dtype=l2h.type, order='F')

        w, h = l2h.width, l2h.height

        rcorr = np.zeros((l2h.N, h, w), dtype=l2h.type)#, order='F').T
        rcorrg = np.zeros((l2h.N, h, w), dtype=l2h.type)#, order='F').T
        aot550pix = np.zeros((w, h), dtype=l2h.type, order='F').T
        betapix = np.zeros((w, h), dtype=l2h.type, order='F').T
        brdfpix = np.zeros((w, h), dtype=l2h.type, order='F').T

        l2h.l2_product.getBand('SZA').writePixels(0, 0, w, h, l2h.sza)
        l2h.l2_product.getBand('VZA').writePixels(0, 0, w, h, np.array(l2h.vza[1]))
        l2h.l2_product.getBand('AZI').writePixels(0, 0, w, h, np.array(l2h.razi[1]))

        ######################################
        #      Add terrain attributes
        ######################################
        if dem:
            sza_mean, sazi_mean = np.nanmean(l2h.sza),np.nanmean(l2h.sazi)
            l2h.slope, l2h.shade = _utils.get_dem_attributes(l2h.elevation, sza=sza_mean, sun_azi=sazi_mean)
            # add elevation band
            l2h.l2_product.getBand('elevation').writePixels(0, 0, w, h, l2h.elevation)
            l2h.l2_product.getBand('slope').writePixels(0, 0, w, h, l2h.slope)
            l2h.l2_product.getBand('shade').writePixels(0, 0, w, h, l2h.shade)

        ######################################
        #      MAIN LOOP
        ######################################
        logging.info('processing ' + file + '...')

        ######################################
        #      First step: AOT adjustment
        ######################################
        # TODO slice/reshape.. to use the two SWIR bands only
        # xblock, yblock = 25, 25
        # TODO aot estimation on water pixels, check better spatial resolution

        ######################################
        #      Second step: Atmosphere and surface correction
        ######################################
        # TODO put chunck size in config yaml file
        for iy in range(0, w, yblock):
            logging.info('process row ' + str(iy) + ' / ' + str(h))
            yc = iy + yblock
            if yc > w:
                yc = w
            for ix in range(0, h, xblock):
                #print('process col ' + str(ix) + ' / ' + str(w))
                xc = ix + xblock
                if xc > h:
                    xc = h
                sza = l2h.sza[ix:xc, iy:yc]
                xshape, yshape = sza.shape

                if (xshape == 0) or (yshape == 0):
                    continue

                razi = l2h.razi[:, ix:xc, iy:yc]
                vza = l2h.vza[:, ix:xc, iy:yc]
                muv = l2h.muv[:, ix:xc, iy:yc]
                mu0 = l2h.mu0[ix:xc, iy:yc]
                mask = l2h.mask[ix:xc, iy:yc]
                flags = l2h.flags[ix:xc, iy:yc]
                band_rad = l2h.band_rad[:, ix:xc, iy:yc]

                maskpixels = maskpixels_
                if allpixels:
                    maskpixels = maskpixels * 0
                elif waterdetect_only:
                    #print('watermask', l2h.watermask[ix:xc, iy:yc].shape)
                    maskpixels[l2h.watermask[ix:xc, iy:yc] == 1] = 0
                else:
                    maskpixels = mask

                if dem:
                    elev = l2h.elevation[ix:xc, iy:yc]
                    pressure = acutils.misc.get_pressure(elev, l2h.pressure_msl)
                    pressure_corr = pressure / l2h.pressure_ref
                else:
                    pressure_corr = np.full((xshape, yshape),l2h.pressure / l2h.pressure_ref)
                pressure_corr = np.array(pressure_corr, dtype=l2h.type, order='F')

                # ---------
                # if maja L2A image provided, use AOT_MAJA product
                # AOT_maja seems to be largely overestimated, before further analyses: force usage of CAMS instead
                # if maja:
                #     aot550guess = l2h.aot_maja[i]
                #     # aot550guess[aot550guess < 0.01] = 0.01
                # else:
                #     aot550guess = aot550rast[i]

                aot550guess = np.array(aot550rast[ix:xc, iy:yc], dtype=l2h.type, order='F')
                fcoef = np.array(fcoefrast[ix:xc, iy:yc], dtype=l2h.type, order='F')
                aot_tot = np.array(aotrast[:, ix:xc, iy:yc], dtype=l2h.type, order='F')
                if l2h.aerosol == 'cds_forecast':
                    aot_sca = np.array(aotscarast[:, ix:xc, iy:yc], dtype=l2h.type, order='F')
                else:
                    aot_sca = aot_tot

                for iband in range(l2h.N):
                    # correct for gaseous absorption
                    tg = smac.compute_gas_trans(iband, l2h.pressure_msl, mu0.reshape(-1),
                                                muv[iband].reshape(-1)).reshape(xshape, yshape)
                    band_rad[iband] = band_rad[iband] / tg



                p = grs_solver.grs.main_algo(xshape, yshape, *lutf.refl.shape,
                                         aotlut, sza_, azi_, vza_,
                                         lutf.refl, lutc.refl, lutf.Cext, lutc.Cext,
                                         vza, sza, razi, band_rad, maskpixels,
                                         l2h.wl, pressure_corr, l2h.sensordata.rg, l2h.solar_irr, l2h.rot,
                                         aot_tot, aot_sca, aot550guess, fcoef,
                                         l2h.nodata, l2h.rrs)


                rcorr[:, ix:xc, iy:yc] = p[0]  # np.transpose(p[0],(0,2,1))
                rcorrg[:, ix:xc, iy:yc] = p[1]  # np.transpose(p[1],(0,2,1))
                aot550pix[ix:xc, iy:yc] = p[2]
                brdfpix[ix:xc, iy:yc] = p[3]

                # TODO improve checksum scheme
                l2h.checksum('row:col, ' + str(ix) + ':' + str(iy))

        rcorr[rcorr == l2h.nodata] = np.nan
        rcorrg[rcorrg == l2h.nodata] = np.nan

        ndwi_corr = np.array((rcorrg[l2h.sensordata.NDWI_vis] - rcorrg[l2h.sensordata.NDWI_nir]) / \
                             (rcorrg[l2h.sensordata.NDWI_vis] + rcorrg[l2h.sensordata.NDWI_nir]))
        # set flags
        l2h.flags = l2h.flags + \
                ((l2h.mask == 1) +
                 (np.array((rcorr[1] < -0.01) | (rcorr[2] < -0.01)) << 1) +
                 ((l2h.mask == 2) << 2) +
                 (((ndwi_corr < l2h.sensordata.NDWI_threshold[0]) | (
                         ndwi_corr > l2h.sensordata.NDWI_threshold[1])) << 3) +
                 ((rcorrg[l2h.sensordata.high_nir[0]] > l2h.sensordata.high_nir[1]) << 4)
                 )
        #logging.info(w, h, l2h.flags.astype(np.uint32).shape,brdfpix.shape,aot550pix.shape,rcorr.shape)
        l2h.l2_product.getBand('flags').writePixels(0, 0, w, h, np.array(l2h.flags.astype(np.uint32)))
        l2h.l2_product.getBand('BRDFg').writePixels(0, 0, w, h, brdfpix)
        l2h.l2_product.getBand("aot550").writePixels(0, 0, w, h, aot550pix)
        for iband in range(l2h.N):
            l2h.l2_product.getBand(l2h.output + '_' + l2h.band_names[iband]). \
                writePixels(0, 0, w, h, rcorr[iband])
            l2h.l2_product.getBand(l2h.output + '_g_' + l2h.band_names[iband]). \
                writePixels(0, 0, w, h, rcorrg[iband])

        logging.info("finishing writing product...")
        l2h.finalize_product()

        logging.info("If no error occured, product is available here : "+outfile)
        if unzip:
            # remove unzipped files (Sentinel files)
            shutil.rmtree(file, ignore_errors=True)

        if untar:
            # remove untared files (Landsat files)
            shutil.rmtree(tmp_dir, ignore_errors=True)
