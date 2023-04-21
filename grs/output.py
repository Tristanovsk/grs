
class l2a():
    def __init__(self,prod):
        self.prod =prod

    def create_product(self, maja=None, waterdetect=None):
        '''
        Create output product dimensions, variables, attributes, flags....
        :return:
        '''

        # allocate flags param
        self.maja = maja
        self.waterdetect = waterdetect

        product = self.product
        ac_product = Product('L2grs', 'L2grs', self.width, self.height)
        writer = ProductIO.getProductWriter('NetCDF4-BEAM')  #
        # writer = ProductIO.getProductWriter('BEAM-DIMAP')
        self.outfile_ext = self.outfile + '.nc'  # dim'
        ac_product.setProductWriter(writer)

        ProductUtils.copyGeoCoding(product, ac_product)
        ProductUtils.copyMetadata(product, ac_product)
        ac_product.setStartTime(product.getStartTime())
        ac_product.setEndTime(product.getEndTime())

        # set metadata: ancillary data used for processing
        meta = jpy.get_type('org.esa.snap.core.datamodel.MetadataElement')
        att = jpy.get_type('org.esa.snap.core.datamodel.MetadataAttribute')
        # att(name=string,type=int), type: 41L->ascii; 12L->int32;

        meta = meta('L2')
        meta.setName('Ancillary Data')
        att_ = att('sensor', ProductData.TYPE_ASCII)
        att_.setDataElems(self.sensor)
        meta.addAttribute(att_)
        att_ = att('processor', ProductData.TYPE_ASCII)
        att_.setDataElems(self.processor)
        meta.addAttribute(att_)
        att_ = att('wkt', ProductData.TYPE_ASCII)
        att_.setDataElems(self.wkt)
        meta.addAttribute(att_)
        # TODO pass int or float argument in proper format instead converting into string
        att_ = att('aerosol_data', ProductData.TYPE_ASCII)
        att_.setDataElems(self.aerosol)
        meta.addAttribute(att_)
        att_ = att('cams_file', ProductData.TYPE_ASCII)
        att_.setDataElems(self.ancillary)
        meta.addAttribute(att_)
        att_ = att('aeronet_file', ProductData.TYPE_ASCII)
        att_.setDataElems(self.aeronetfile)
        meta.addAttribute(att_)
        att_ = att('aot_spectral', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aot))
        meta.addAttribute(att_)
        att_ = att('aot550', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aot550))
        meta.addAttribute(att_)
        att_ = att('pressure_sea_level', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.pressure_msl))
        att_.setUnit('hPa')
        meta.addAttribute(att_)
        att_ = att('ozone', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.o3du))
        att_.setUnit('DU')
        meta.addAttribute(att_)
        att_ = att('water_vapor', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.h2o))
        att_.setUnit('g cm-2')
        meta.addAttribute(att_)
        att_ = att('no2', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.tcno2))
        meta.addAttribute(att_)
        att_ = att('wind_module', ProductData.TYPE_ASCII)
        att_.setDataElems(str('...'))
        meta.addAttribute(att_)
        att_ = att('temperature_2m', ProductData.TYPE_ASCII)
        att_.setDataElems(str(self.aux.t2m))
        meta.addAttribute(att_)

        ac_product.getMetadataRoot().addElement(meta)

        # -------------------------------------
        # set masks / flags
        # -------------------------------------

        flags = ac_product.addBand('flags', ProductData.TYPE_UINT32)
        flags.setDescription('Flags for aquatic color purposes')
        # vflags = ac_product.addBand('valid', ProductData.TYPE_UINT8)
        # vflags.setDescription('used to set valid data pixels')
        expr_valid_pixel = 'mask_nodata == 0 && mask_negative == 0'

        # Also for each flag a layer should be created
        Color = jpy.get_type('java.awt.Color')
        colors = [Color.BLUE, Color.YELLOW, Color.RED, Color.PINK, Color.MAGENTA, Color.GREEN, Color.GRAY] * 10
        coding = FlagCoding('flags')
        f = [None] * 32
        mask_id = 0
        f[mask_id] = coding.addFlag("nodata", 2 ** mask_id, "nodata in input image ")
        mask_id += 1
        f[mask_id] = coding.addFlag("negative", 2 ** mask_id, "negative values in visible ")
        mask_id += 1
        f[mask_id] = coding.addFlag("ndwi", 2 ** mask_id, "based on ndwi vis nir TOA based on bands " +
                                    self.band_names[self.sensordata.NDWI_vis] + " and " +
                                    self.band_names[self.sensordata.NDWI_nir] +
                                    " for range " + str(self.sensordata.NDWI_threshold))
        mask_id += 1
        f[mask_id] = coding.addFlag("ndwi_swir", 2 ** mask_id, "based on ndwi nir swir TOA based on bands" +
                                    self.band_names[self.sensordata.NDWI_swir_nir] + " and " +
                                    self.band_names[self.sensordata.NDWI_swir_swir] +
                                    " for range " + str(self.sensordata.NDWI_swir_threshold))
        mask_id += 1
        f[mask_id] = coding.addFlag("high_nir", 2 ** mask_id,
                                    "high radiance in the nir band (e.g., cloud, snow); condition Rrs_g at " +
                                    self.band_names[self.sensordata.high_nir[0]] + " greater than " + str(
                                        self.sensordata.high_nir[1]))
        mask_id += 1
        f[mask_id] = coding.addFlag("hicld", 2 ** mask_id,
                                    "high cloud as observed from cirrus band; condition Rtoa at band " +
                                    self.sensordata.cirrus[0] + " greater than " + str(self.sensordata.cirrus[1]))
        if self.sensordata.O2band:
            O2info1 = self.sensordata.O2band[0] + " greater than " + str(self.sensordata.O2band[1])
            O2info2 = self.sensordata.O2band[0] + " greater than " + str(self.sensordata.O2band[2])
        else:
            O2info1 = 'non applicable'
            O2info2 = O2info1

        mask_id += 1
        f[mask_id] = coding.addFlag("moderate_cloud_risk_O2band", 2 ** mask_id,
                                    "moderate risk of bright cloud as observed from O2 band; condition Rtoa at band " +
                                    O2info1)
        mask_id += 1
        f[mask_id] = coding.addFlag("high_cloud_risk_O2band", 2 ** mask_id,
                                    "high risk of bright cloud as observed from O2 band; condition Rtoa at band " +
                                    O2info2)
        mask_id += 1
        f[mask_id] = coding.addFlag("L1_opaque_clouds", 2 ** mask_id, " flag from L1 image ")
        mask_id += 1
        f[mask_id] = coding.addFlag("L1_cirrus", 2 ** mask_id, "cirrus cloud flag from L1 image ")
        mask_id += 1
        f[mask_id] = coding.addFlag("L1_shadow", 2 ** mask_id, "cloud-shadow flag from L1 image ")

        for i_f in range(mask_id):
            ac_product.addMask('mask_' + f[i_f].getName(), 'flags.' + f[i_f].getName(),
                               f[i_f].getDescription(), colors[i_f], 0.3)

        # -------------------
        # for Sentinel 2
        # if MAJA L2A / WaterDetect provided, load respective flags
        # WARNING: mask_id must remain smaller than 32 (binary coding)
        mask_id += 1
        additional_f = []
        # if self.maja:

        for i, mask in enumerate(self.maja_masks):
            logging.info('Mask binary ' + str(mask_id))
            additional_f.append(coding.addFlag(mask, 2 ** mask_id, 'Mask ' + mask + ' imported from MAJA chain'))
            mask_id += 1

        if self.waterdetect:
            logging.info('waterdetect mask ' + mask_id)
            additional_f.append(
                coding.addFlag('WaterDetect', 2 ** mask_id, 'Water mask imported from WaterDetect processing'))

        ac_product.getFlagCodingGroup().add(coding)
        flags.setSampleCoding(coding)

        # -------------------
        # for Sentinel 2
        # if MAJA L2A / WaterDetect provided, load respective flags
        # WARNING: mask_id must remain smaller than 32 (binary coding)

        mask_id = 0
        # to get fixed format of output image put the MAJA flags even if not available
        # if self.maja:

        for i, mask in enumerate(self.maja_masks):
            f = additional_f[mask_id]
            ac_product.addMask('' + f.getName(), 'flags.' + f.getName(),
                               f.getDescription(), colors[mask_id], 0.3)
            mask_id += 1

        if self.waterdetect:
            f = additional_f[mask_id]
            ac_product.addMask('mask_' + f.getName(), 'flags.' + f.getName(),
                               f.getDescription(), colors[mask_id], 0.3)

        # set data
        # Water-leaving radiance + sunglint
        for iband in range(self.N):
            bname = self.output + '_g_' + self.band_names[iband]
            acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
            acband.setSpectralWavelength(self.wl[iband])
            acband.setSpectralBandwidth(self.B[iband].getSpectralBandwidth())
            acband.setModified(True)
            acband.setNoDataValue(np.nan)
            acband.setNoDataValueUsed(True)
            acband.setValidPixelExpression(expr_valid_pixel)
            if self.output == 'Lwn':
                ac_product.getBand(bname).setDescription(
                    'Water-leaving plus sunglint normalized radiance (Lwn + Lg) in mW cm-2 sr-1 μm-1 at ' +
                    self.band_names[iband])
            else:
                ac_product.getBand(bname).setDescription(
                    'Water-leaving plus sunglint remote sensing reflectance (Rrs + Lg/F0) in sr-1 at ' +
                    self.band_names[iband])

        # Water-leaving radiance
        for iband in range(self.N):
            bname = self.output + '_' + self.band_names[iband]
            acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
            acband.setSpectralWavelength(self.wl[iband])
            acband.setSpectralBandwidth(self.B[iband].getSpectralBandwidth())
            acband.setModified(True)
            acband.setNoDataValue(np.nan)
            acband.setNoDataValueUsed(True)
            acband.setValidPixelExpression(expr_valid_pixel)
            if self.output == 'Lwn':
                ac_product.getBand(bname).setDescription(
                    'Normalized water-leaving radiance in mW cm-2 sr-1 μm-1 at ' + self.band_names[iband])
            else:
                ac_product.getBand(bname).setDescription(
                    'Remote sensing reflectance in sr-1 at ' + self.band_names[iband])
        # Sunglint reflection factor
        # for iband in range(self.N):
        bname = 'BRDFg'  # + self.band_names[iband]
        acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
        # acband.setSpectralWavelength(self.wl[iband])
        # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        acband.setValidPixelExpression(expr_valid_pixel + ' && ' + bname + ' >= 0')
        ac_product.getBand(bname).setDescription('Glint reflection factor (BRDF) ')  # + self.band_names[iband])

        # estimated aerosol optical thickness at 550 nm
        bname = 'aot550'  # + self.band_names[iband]
        acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
        # acband.setSpectralWavelength(self.wl[iband])
        # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        acband.setValidPixelExpression(expr_valid_pixel + ' && ' + bname + ' >= 0')
        ac_product.getBand(bname).setDescription('aerosol optical thickness at 550 nm ')  # + self.band_names[iband])

        # if maja option is on, copy aerosol optical thickness product from MAJA L2A image
        # if self.maja:
        bname = 'aot_maja'  # + self.band_names[iband]
        acband = ac_product.addBand(bname, ProductData.TYPE_FLOAT32)
        # acband.setSpectralWavelength(self.wl[iband])
        # acband.setSpectralBandwidth(self.b[iband].getSpectralBandwidth())
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        acband.setValidPixelExpression(expr_valid_pixel + ' && ' + bname + ' >= 0')
        ac_product.getBand(bname).setDescription(
            'AOT product from MAJA processing (L2A)')  # + self.band_names[iband])

        # Viewing geometry
        acband = ac_product.addBand('SZA', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('SZA').setDescription('Solar zenith angle in deg.')

        acband = ac_product.addBand('VZA', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('VZA').setDescription('Mean viewing zenith angle in deg.')

        acband = ac_product.addBand('AZI', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('AZI').setDescription('Mean relative azimuth angle in deg.')

        acband = ac_product.addBand('ndwi', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('ndwi').setDescription("based on ndwi vis nir TOA based on bands " +
                                                  self.band_names[self.sensordata.NDWI_vis] + " and " +
                                                  self.band_names[self.sensordata.NDWI_nir])

        acband = ac_product.addBand('ndwi_swir', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('ndwi_swir').setDescription("based on ndwi nir swir TOA based on bands " +
                                                       self.band_names[self.sensordata.NDWI_swir_nir] + " and " +
                                                       self.band_names[self.sensordata.NDWI_swir_swir])

        acband = ac_product.addBand('elevation', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        # TODO update the name of DEM data set used in function of the input (e.g., SRTM, GETASSE)
        ac_product.getBand('elevation').setDescription('elevation from SRTM 3Sec in meters')

        acband = ac_product.addBand('shade', ProductData.TYPE_FLOAT32)
        acband.setModified(True)
        acband.setNoDataValue(np.nan)
        acband.setNoDataValueUsed(True)
        ac_product.getBand('shade').setDescription('terrain shades computed from DEM and sun geometry')

        ac_product.setAutoGrouping(self.output + ':' + self.output + '_g_')

        logging.info("Creating output folder "+self.outfile_ext)
        ac_product.writeHeader(String(self.outfile_ext))
        # next line needed since snap 'writeHeader' force the extension to be consistent with data type (e.g., .tif for GeoTIFF)
        os.rename(self.outfile_ext, self.outfile_ext + '.incomplete')

        self.l2_product = ac_product

    def checksum(self, info):
        '''
        Save info on the current processing stage
        :param info:
        :return:
        '''
        # TODO improve checksum scheme
        with open(self.outfile + '.checksum', "w") as f:
            f.write(info)

    def deleteCache(self):
        S2CacheUtils = jpy.get_type('org.esa.s2tbx.dataio.cache.S2CacheUtils')
        S2CacheUtils.deleteCache()

    def finalize_product(self):
        '''remove checksum file
        remove extension ".incomplete" from output file name
        convert into netcdf (compressed) from gpt and ncdump/nco tool'''
        # jpy.destroy_jvm()
        self.product.dispose()
        self.l2_product.dispose()

        os.remove(self.outfile + '.checksum')
        name = self.outfile_ext + '.incomplete'
        # final_name = os.path.splitext(name)[0]
        os.rename(name, self.outfile_ext)  # final_name)

    def print_info(self):
        ''' logging.info info, can be used to check if object is complete'''
        logging.info('Product: %s, %d x %d pixels, %s' % (self.name, self.width, self.height, self.description))
        logging.info('Bands:   %s' % (list(self.band_names)))
        for i in range(len(self.wl)):
            logging.info('Band ' + str(i) + ", " + str(self.B[i].getSpectralWavelength()) + ' centered on ' + str(
                self.wl[i]) + 'nm loaded')
