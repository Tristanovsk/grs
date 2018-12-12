from pathlib import Path
from esasnappy import ProductData, ProductIO

from . import acutils
from . import auxdata
from . import utils
from .anglegen import *
from .fortran import main_algo as f


class process:
    def __init__(self):
        pass

    def execute(self, file, sensor, wkt, altitude=0, aerosol='cams', noclobber=True, outfile=None,
                gdm=None, aeronet_file=None, aot550=0.1, angstrom=1, resolution=None, indband=None):
        '''

        :param file:
        :param sensor:
        :param wkt:
        :param output:
        :param gdm:
        :param altitude:
        :param noclobber:
        :param aeronet_file:
        :param resolution:
        :param bands_selection:
        :return:
        '''

        ##################################
        # File naming convention
        ##################################
        if outfile == None:
            if 'S2' in sensor:
                outfile = file.replace('L1C', 'L2h')
                outfile = outfile.replace('.SAFE', '').rstrip('/')
            elif 'LANDSAT' in sensor:
                outfile = file.replace('L1TP', 'L2h')
                outfile = outfile.replace('.txt', '').rstrip('/')
            else:
                print('Not recognized sensor, please try again!')
                sys.exit()

        basename = os.path.basename(outfile)

        if os.path.isfile(outfile + ".dim") & os.path.isdir(outfile + ".data") & noclobber:
            print('File ' + outfile + ' already processed; skip!')
            sys.exit()

        ##################################
        # Get sensor auxiliary data
        ##################################
        sensordata = auxdata.sensordata(sensor)
        if resolution == None:
            resolution = sensordata.resolution
        indband = sensordata.indband

        ##################################
        # Read L1C product
        ##################################
        print("Reading...")
        print(file)
        product = ProductIO.readProduct(file)

        ##################################
        # Generate l2h object
        ##################################
        _utils = utils.utils()
        l2h = utils.info(product, sensordata, aerosol)

        ##################################
        # Set bands to be processed including NIR and SWIR
        ##################################
        l2h.headerfile = file
        l2h.band_names = l2h.sensordata.band_names

        ##################################
        # GET METADATA
        ##################################
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
        # GENERATE BAND ANGLES
        ##################################
        if 'S2' in sensor:
            # TODO remove next line when the S2resampler is fixed up by ESA (hopefully soon, 05/10/2018)
            # angle_generator().sentinel2(l2h)
            pass
        elif 'LANDSAT_8' in sensor:
            angle_generator().landsat(l2h)
        else:
            angle_generator().landsat_tm(l2h)

        ##################################
        # RESAMPLE TO A UNIQUE RESOLUTION
        ##################################
        if 'S2' in sensor:
            l2h.product = _utils.s2_resampler(l2h.product, resolution=resolution)
            # l2h.product = _utils.generic_resampler(l2h.product, resolution=resolution, method='Nearest')
        else:
            l2h.product = _utils.resampler(l2h.product, resolution=resolution)  # , upmethod='Nearest')

        ##################################
        # SUBSET TO AREA OF INTEREST
        ##################################

        l2h.product = _utils.get_subset(l2h.product, wkt)
        l2h.get_product_info()
        l2h.set_outfile(outfile)
        l2h.wkt = _utils.get_extent(l2h.product)
        l2h.crs = str(l2h.product.getBand(l2h.band_names[0]).getGeoCoding().getImageCRS())

        ##################################
        # GET IMAGE AND RASTER PROPERTIES
        ##################################

        l2h.get_bands(l2h.band_names)
        l2h.print_info()

        ##################################
        # SET NEW BAND FOR MASKING
        ##################################
        l2h.product.addBand('ndwi', ProductData.TYPE_FLOAT32)
        l2h.ndwi_band = l2h.product.getBand('ndwi')
        l2h.ndwi_band.ensureRasterData()
        l2h.ndwi_band.loadRasterData()

        ##################################
        # GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
        ##################################
        l2h.aux = auxdata.cams()
        target = Path(os.path.join(l2h.cams_folder, l2h.date.strftime('%Y-%m') + '_month_' + l2h.aerosol + '.nc'))
        l2h.aux.load_cams_data(target, l2h.date, data_type=l2h.aerosol)
        l2h.aux.get_cams_ancillary(target, l2h.date, l2h.wkt)

        ## uncomment this part to use ecmwf files provided in the .SAFE format
        # if 'S2' in sensor:
        #     l2h.aux.get_tile_dir(file)
        #     l2h.aux.get_aux_dir()
        #     l2h.aux.get_ecmwf_data()

        # get pressure at the scene altitude
        l2h.pressure = acutils.misc.get_pressure(altitude, l2h.aux.msl)
        l2h.aux.pressure = l2h.pressure

        ##################################
        # GET ANCILLARY DATA (AEROSOL)
        ##################################
        aero = acutils.aerosol()

        # AERONET data
        if (l2h.aerosol == 'aeronet'):
            l2h.set_aeronetfile(aeronet_file)
            try:
                l2h.aux.Aeronet.import_aeronet_data(aero, l2h.aeronetfile, l2h.date)
            except:
                print('Error: No aeronet data in the +/- 2day time window.')
                sys.exit()

            l2h.aux.aot_wl = aero.wavelengths
            l2h.aux.aot = aero.aot
            l2h.aux.aot550 = aero.aot550

        # CAMS (ECMWF)
        elif (l2h.aerosol == 'cams_forecast') | (l2h.aerosol == 'cams_reanalysis'):
            # daily file
            # target=Path(os.path.join(l2h.ecmwf_root,l2h.date.strftime('%Y-%m-%d')+'_cams_aero.nc'))
            # cams=aux.get_cams_aerosol(target,l2h.date.strftime('%Y-%m-%d'),l2h.wkt,l2h.crs)
            # monthly file

            l2h.aux.get_cams_aerosol(target, l2h.date, l2h.wkt)

        elif (l2h.aerosol == 'user_model'):
            l2h.aot550 = aot550
            l2h.angstrom = angstrom
            l2h.aux.aot = l2h.aot550 * (np.array(l2h.wl) / 550) ** (-l2h.angstrom)
            l2h.aux.aot_wl = l2h.wl

        else:
            sys.exit("No aerosol data provided, try again.")

        # set spectral aot for satellie bands
        aero.fit_spectral_aot(l2h.aux.aot_wl, l2h.aux.aot)
        l2h.aot = aero.get_spectral_aot(np.array(l2h.wl))
        l2h.aot550 = l2h.aux.aot550

        #####################################
        # LOAD LUT FOR ATMOSPHERIC CORRECTION
        #####################################
        # load lut
        print('loading lut...', l2h.lutfine)
        lutf = acutils.lut(l2h.band_names)
        lutc = acutils.lut(l2h.band_names)
        lutf.load_lut(l2h.lutfine, indband)
        lutc.load_lut(l2h.lutcoarse, indband)

        # normalization of Cext to get spectral dependence of fine and coarse modes
        nCext_f = lutf.Cext / lutf.Cext550
        nCext_c = lutc.Cext / lutc.Cext550
        aero.fit_aero(nCext_f, nCext_c, l2h.aot / l2h.aot550)
        l2h.fcoef = aero.fcoef
        # rlut = [x + y for x, y in zip([fcoef * x for x in rlut_f], [(1. - fcoef) * x for x in rlut_c])]

        # acutils.plot_lut(grid_lut[3], grid_lut[2], rlut_c[9][2, 13, :, :])
        # acutils.plot_lut(grid_lut[3], grid_lut[2], rlut_f[9][2, 13, :, :])
        # acutils.plot_lut(grid_lut[3], grid_lut[2], rlut[9][2, 13, :, :])

        # correction for pressure
        # scaling of the lut values and rayleigh optical thickness
        rlut_f = [x * l2h.pressure / l2h.pressure_ref for x in lutf.refl]
        rlut_c = [x * l2h.pressure / l2h.pressure_ref for x in lutc.refl]
        l2h.rot = [x * l2h.pressure / l2h.pressure_ref for x in l2h.sensordata.rot]

        ####################################
        #     Set SMAC patrameters for
        #    absorbing gases correction
        ####################################
        smac = acutils.smac(l2h.sensordata.smac_bands, l2h.sensordata.smac_dir)
        smac.set_gas_param()
        smac.set_values(o3du=l2h.aux.o3du, h2o=l2h.aux.h2o)
        smac.set_standard_values(l2h.pressure)
        l2h.aux.no2 = smac.uno2

        ######################################
        #      Create output l2 product
        #          'l2_product'
        ######################################
        l2h.create_product()
        # reshaping for fortran binding
        aotlut = np.array(lutf.aot, dtype=l2h.type, order='F')
        vzalut = np.array(lutf.vza, dtype=l2h.type, order='F')
        szalut = np.array(lutf.sza, dtype=l2h.type, order='F')
        razilut = np.array(lutf.azi, dtype=l2h.type, order='F')
        rlut_f = np.array(rlut_f, dtype=l2h.type, order='F')
        rlut_c = np.array(rlut_c, dtype=l2h.type, order='F')
        grid_lut = (szalut, razilut, vzalut)

        ######################################
        #      MAIN LOOP
        ######################################
        # arrays allocation
        aot550pix = np.zeros(l2h.width, dtype=l2h.type)
        rtoaf = np.zeros((lutf.aot.__len__(), l2h.N, l2h.width), dtype=l2h.type, order='F')
        rtoac = np.zeros((lutc.aot.__len__(), l2h.N, l2h.width), dtype=l2h.type, order='F')

        # set aot by hand
        aot550pix.fill(l2h.aot550)

        for i in range(l2h.height):
            print('process row ' + str(i))
            # LOAD PIXELS DATA FOR ROW #i
            l2h.load_data(i)

            validx = (l2h.mask == 0)
            # l2h.sza[~validx]=np.nan
            # l2h.razi[~validx]=np.nan

            sza = l2h.sza[validx]
            razi = l2h.razi[validx]
            vza = l2h.vza[validx]

            for iband in range(l2h.N):
                # preparing lut data
                grid_pix = list(zip(sza, razi[..., iband], vza[..., iband]))

                for iaot in range(aotlut.__len__()):
                    rtoaf[iaot, iband, validx] = lutf.interp_lut(grid_lut, rlut_f[iband][iaot, ...], grid_pix)
                    rtoac[iaot, iband, validx] = lutc.interp_lut(grid_lut, rlut_c[iband][iaot, ...], grid_pix)

                # correct for gaseous absorption
                tg = smac.compute_gas_trans(iband, l2h.pressure, l2h.mu0, l2h.muv[iband])
                l2h.rs2[:, iband] = l2h.rs2[:, iband] / tg

            rcorr, rcorrg = f.main_algo(l2h.width, l2h.N, aotlut.__len__(),
                                        l2h.vza, l2h.sza, l2h.razi, l2h.rs2, l2h.mask, l2h.wl,
                                        aotlut, rtoaf, rtoac, lutf.Cext, lutc.Cext,
                                        l2h.sensordata.rg, l2h.solar_irr, l2h.rot,
                                        l2h.aot, aot550pix, l2h.fcoef, l2h.nodata)

            # reshape for snap modules
            rcorr = np.ma.array(rcorr.T, mask=rcorr.T == l2h.nodata, fill_value=np.nan)  # .tolist()
            rcorrg = np.ma.array(rcorrg.T, mask=rcorrg.T == l2h.nodata, fill_value=np.nan)  # .tolist()

            ndwi_corr = np.array((rcorrg[l2h.sensordata.NDWI_vis] - rcorrg[l2h.sensordata.NDWI_nir]) / \
                                 (rcorrg[l2h.sensordata.NDWI_vis] + rcorrg[l2h.sensordata.NDWI_nir]))
            # set flags
            l2h.flags = ((l2h.mask == 1) +
                         (np.array((rcorr[1] < -0.01) | (rcorr[2] < -0.01)) << 1) +
                         ((l2h.mask == 2) << 2) +
                         (((ndwi_corr < l2h.sensordata.NDWI_threshold[0]) | (
                                 ndwi_corr > l2h.sensordata.NDWI_threshold[1])) << 3) +
                         ((rcorrg[l2h.sensordata.NDWI_nir] > 5.) << 4)
                         )

            for iband in range(l2h.N):
                l2h.l2_product.getBand("Lwn_" + l2h.band_names[iband]). \
                    writePixels(0, i, l2h.width, 1, rcorr[iband])
                l2h.l2_product.getBand("Lwn_g_" + l2h.band_names[iband]). \
                    writePixels(0, i, l2h.width, 1, rcorrg[iband])

            l2h.l2_product.getBand('flags').writePixels(0, i, l2h.width, 1, np.array(l2h.flags, dtype=np.uint32))
            l2h.l2_product.getBand("BRDFg").writePixels(0, i, l2h.width, 1, l2h.ndwi)
            l2h.l2_product.getBand("SZA").writePixels(0, i, l2h.width, 1, l2h.sza)
            l2h.l2_product.getBand("VZA").writePixels(0, i, l2h.width, 1, np.array(l2h.vza[:, 1]))
            l2h.l2_product.getBand("AZI").writePixels(0, i, l2h.width, 1, np.array(l2h.razi[:, 1]))
