import os

import importlib_resources
import yaml

import numpy as np
import xarray as xr

# keep attributes through operatyion on xarray objects
xr.set_options(keep_attrs=True)
import rioxarray as rio
import logging
import gc

from multiprocessing import Pool  # Process pool
from multiprocessing import sharedctypes
import itertools

import GRSdriver

from . import Product, acutils, AuxData, CamsProduct, L2aProduct, Masking, Rasterization, Kernel

opj = os.path.join

configfile = importlib_resources.files(__package__) / 'config.yml'
with open(configfile, 'r') as file:
    config = yaml.safe_load(file)

GRSDATA = config['path']['grsdata']
TOALUT = config['path']['toa_lut']
TRANSLUT = config['path']['trans_lut']
CAMS_PATH = config['path']['trans_lut']
NCPU = config['processor']['ncpu']
NETCDF_ENGINE = config['processor']['netcdf_engine']


class Process:
    '''
    Main GRS class.

    '''

    def __init__(self):
        self.lut_file = opj(GRSDATA, TOALUT)
        self.trans_lut_file = opj(GRSDATA, TRANSLUT)
        self.cams_dir = CAMS_PATH
        self.Nproc = NCPU
        self.pressure_ref = 101500.
        self.flags_tokeep = [3]
        self.flags_tomask = [0, 1, 10, 13, 14, 18]
        self.successful = False

    def execute(self, l1c_prod,
                ofile='',
                cams_file=None,
                surfwater_file=None,
                dem_file=None,
                resolution=20,
                megapix_size=4,
                chunk = 512,
                aot_megapix_size=40,
                aot_chunk = 25,
                scale_aot=1,
                opac_model=None,
                allpixels=False,
                snap_compliant=False,
                provide_kernel=False
                ):

        '''
        Main program calling all GRS steps

        :param l1c_prod: xarray L1C object or L1C input file (path) to be processed
        :param ofile: Absolute path of the output file
        :param cams_file: Absolute path for root directory of CAMS data
        :param surfwater_file: Absolute path the surfwater file (.tif)
        :param dem_file: Absolute path of the DEM geotiff file
        :param resolution: pixel resolution in meter
        :param scale_aot: scaling factor applied to CAMS aod550 raster
        :param opac_model: If set force OPAC aerosol model for LUT interpolation (taken from CAMS data otherwise)
        :param allpixels: if True process all pixels (no water pixel masking)
        :param snap_compliant: Output format compliant with SNAP software for practical analysis
        :return:

        Examples
        --------

        >>> import grs
        >>> file ='$YOUR_PATH_TO_IMG/S2A_MSIL1C_20181019T102031_N0500_R065_T30PYT_20230815T043850.SAFE'
        >>> tile = file.split('_')[-2][1:]
        >>> dem_file = '$YOUR_PATH_TO_DEM/COP-DEM_GLO-30-DGED_'+tile+'.tif'
        >>> cams_file = '$YOUR_PATH_TO_CAMS/cams_forecast_2018-10.nc'
        >>> process_ = grs.Process()

        In this example, you force the aerosol model to be 'DESE_rh70' for desert dust

        >>> process_.execute(file_nc,
        ...                 cams_file=cams_file,
        ...                 surfwater_file=None,
        ...                 dem_file=dem_file,
        ...                 scale_aot=1,
        ...                 opac_model='DESE_rh70')
        INFO:root:pass netcdf image as grs product object
        INFO:root:get CAMS auxilliary data
        INFO:root:flagging from l1c data
        INFO:root:cloud masking with s2cloudless
        INFO:root:land masking
        INFO:root:cirrus masking
        INFO:root:high swir masking
        INFO:root:loading look-up tables
        INFO:root:compute gaseous transmittance from cams data
        INFO:root:correct for gaseous absorption
        INFO:root:compute spectral index (e.g., NDWI)
        INFO:root:apply water masking
        INFO:root:lut interpolation
        INFO:root:selected aerosol model: DESE_rh70
        INFO:root:scaling aot by: 1
        INFO:root:set final parameters
        INFO:root:compute surface pressure from dem
        INFO:root:run grs process
        INFO:root:success
        INFO:root:construct final product
        INFO:root:construct l2a

        >>> process_.l2a.l2_prod
        <xarray.Dataset>
        Dimensions:      (wl: 11, y: 1818, x: 2523)
        Coordinates:
          * wl           (wl) int64 443 490 560 665 705 740 783 842 865 1610 2190
            time         datetime64[ns] 2018-10-19T10:20:31.024000
          * x            (x) float64 7.299e+05 7.299e+05 7.3e+05 ... 7.803e+05 7.804e+05
          * y            (y) float64 1.3e+06 1.3e+06 1.3e+06 ... 1.264e+06 1.264e+06
            band         int64 1
            spatial_ref  int64 0
        Data variables:
            Rrs          (wl, y, x) float32 nan nan nan nan nan ... nan nan nan nan nan
            BRDFg        (y, x) float32 nan nan nan nan nan nan ... nan nan nan nan nan
            aot550       (y, x) float32 0.187 0.187 0.187 0.187 ... 0.1929 0.1929 0.1929
            vza          (y, x) float32 6.489 6.489 6.483 6.483 ... 2.13 2.125 2.125
            sza          (y, x) float32 nan nan nan nan nan nan ... nan nan nan nan nan
            raa          (y, x) float32 332.5 332.5 332.5 332.5 ... 290.0 289.9 289.9
            flags    (y, x) int64 184 56 184 184 184 184 ... 160 160 160 160 160 160
            dem          (y, x) float32 282.9 282.7 282.5 282.3 ... 237.3 237.3 237.4
            surfwater    (y, x) int8 1 1 1 1 1 1 1 1 1 1 1 1 ... 1 1 1 1 1 1 1 1 1 1 1 1
        Attributes: (12/71)
            long_name:                           CA BLUE GREEN RED VRE_1 VRE_2 VRE_3 ...
            constellation:                       Sentinel-2
            constellation_id:                    S2
            product_path:                        /data/satellite/Sentinel-2/L1C/30PYT...
            product_name:                        S2A_MSIL1C_20181019T102031_N0500_R06...
            product_filename:                    S2A_MSIL1C_20181019T102031_N0500_R06...
            ...                                  ...
            ndwi_threshold:                      0.0
            vis_swir_index_threshold:            0.0
            hcld_threshold:                      0.003
            dirdata:                             /data/grs/grsdata
            abs_gas_file:                        /home/harmel/Dropbox/Dropbox/work/gi...
            water_vapor_transmittance_file:      /home/harmel/Dropbox/Dropbox/work/gi...

        You can either further play with the l2a xarray or save it into netcdf:


        >>> process_.ofile='./name_of_your_output_l2a_netcdf'
        >>> process_.write_output()
        INFO:root:export final product into netcdf
        INFO:root:export into encoded netcdf

        '''

        self.ofile = ofile
        self.snap_compliant = snap_compliant
        self.successful = False

        ##################################
        # Get image data
        ##################################
        if isinstance(l1c_prod, str):
            # get extension
            extension = l1c_prod.split('.')[-1]
            basename = os.path.basename(l1c_prod)
            if extension == 'nc':
                logging.info('pass netcdf image as grs product object')
                prod = Product(xr.open_dataset(l1c_prod, engine=NETCDF_ENGINE))
            elif 'SAFE' in extension:
                logging.info('Open L1C Sentinel 2 image and compute angle parameters')
                global l1c
                l1c = GRSdriver.Sentinel2Driver(l1c_prod, resolution=resolution)
                l1c.load_product()
                logging.info('pass raw image as grs product object')
                prod = Product(l1c.prod)
                # clear memory (TODO make it work!!)
                del l1c
                gc.collect()
            elif ('LC09_L1' in basename) or ('LC08_L1' in basename):
                logging.info('Open L1TP Landsat image')

                l1c = GRSdriver.LandsatDriver(l1c_prod, resolution=resolution)
                l1c.load_mask()
                l1c.load_product()
                logging.info('pass raw image as grs product object')
                prod = Product(l1c.prod)
                # clear memory (TODO make it work!!)
                del l1c
                gc.collect()
            else:
                logging.info('input file format not recognized, stop')
                return
        elif isinstance(l1c_prod, xr.Dataset):
            try:
                prod = Product(l1c_prod)
            except:
                logging.info('input file format not recognized, stop')
                return

        self.prod = prod

        ##################################
        # Set sensor specifications
        ##################################
        # TODO check evolution concerning viewing angles computation for Lansdat, now in monoview mode

        if 'S2' in prod.sensor:
            monoview = False
        else:
            monoview = True
        _R_ = Rasterization(monoview=monoview)

        ##################################
        # GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
        ##################################
        logging.info('get CAMS auxilliary data')
        if cams_file:
            cams = CamsProduct(prod.raster, cams_file=cams_file)
        else:
            tile = prod.raster.attrs['tile']
            cams_dir = os.path.join(self.cams_dir, tile)
            cams = CamsProduct(prod.raster, dir=cams_dir, suffix='_' + tile)
        cams.load(daily_stats=False)

        ##################################
        # Pixel classification
        # Generate the flags raster
        ##################################
        logging.info('flagging from l1c data')

        if surfwater_file:
            logging.info('loading surfwater data file')
            prod.raster['surfwater'] = rio.open_rasterio(surfwater_file
                                                         ).astype(np.uint8
                                                                  ).squeeze().interp(x=prod.x,
                                                                                     y=prod.y,
                                                                                     method='nearest')
            prod.raster.surfwater.name = 'surfwater'
            prod.raster.surfwater.attrs = {
                'description': 'surfwater file not provided as input, all pixels flagged as water (e.g., surfwater=1)'}

        masking_ = Masking(prod.raster)
        prod.raster = masking_.process(output="prod")

        # -- clean up
        prod.raster = prod.raster.drop_vars(["surfwater"])

        #####################################
        # SUBSET RASTER TO KEEP REQUESTED BANDS
        #####################################
        # TODO check if we can remove cirrus and water vapor band from output object
        if prod.bcirrus:
            prod.cirrus = prod.raster.bands.sel(wl=prod.bcirrus, method='nearest')
        if prod.bwv:
            prod.wv = prod.raster.bands.sel(wl=prod.bwv, method='nearest')

        prod.raster = prod.raster.sel(wl=prod.wl_process, method='nearest')
        wl_true = prod.raster.wl_true

        ##################################
        ## ADD ELEVATION AND PRESSURE BAND
        ##################################
        # TODO activate DEM loading to improve pressure computation
        # prod.get_elevation()

        #####################################
        # LOAD LUT FOR ATMOSPHERIC CORRECTION
        #####################################
        logging.info('loading look-up tables')
        trans_lut = xr.open_dataset(self.trans_lut_file, engine=NETCDF_ENGINE)
        trans_lut['wl'] = trans_lut['wl'] * 1000

        aero_lut = xr.open_dataset(self.lut_file, engine=NETCDF_ENGINE)
        aero_lut['wl'] = aero_lut['wl'] * 1000
        aero_lut['aot'] = aero_lut.aot.isel(wind=0).squeeze()


        ####################################
        #    absorbing gases correction
        ####################################
        logging.info('compute gaseous transmittance from cams data')

        gas_trans = acutils.GaseousTransmittance(prod, cams)
        gases = ['co2', 'o2', 'o4', 'ch4', 'no2', 'o3', 'h2o']
        for gas in gases:
            gas_trans.coef_abs_scat[gas] = 1
        Tg_raster = gas_trans.get_gaseous_transmittance(gases=['o3', 'no2'])

        logging.info('correct for gaseous absorption')
        for wl in prod.raster.wl.values:
            prod.raster['bands'].loc[wl] = prod.raster.bands.sel(wl=wl) / Tg_raster.sel(wl=wl).interp(x=prod.raster.x,
                                                                 y=prod.raster.y)
        prod.raster.bands.attrs['gas_absorption_correction'] = True

        ######################################
        # Water mask
        ######################################
        # TODO remove ndwi export / replace this part with flags masking instead
        logging.info('compute spectral index (e.g., NDWI)')

        vis = prod.raster.bands.sel(wl=prod.bvis, method='nearest')
        nir = prod.raster.bands.sel(wl=prod.bnir, method='nearest')
        swir = prod.raster.bands.sel(wl=prod.bswir, method='nearest')
        swir2 = prod.raster.bands.sel(wl=prod.bswir2, method='nearest')

        ndwi = (vis - nir) / (vis + nir)
        ndwi_swir = (vis - swir) / (vis + swir)

        prod.raster['ndwi'] = ndwi
        prod.raster.ndwi.attrs = {
            'description': 'Normalized difference spectral index between bands at ' + str(prod.bvis) + ' and ' + str(
                prod.bnir) + ' nm', 'units': '-'}
        prod.raster['ndwi_swir'] = ndwi_swir
        prod.raster.ndwi_swir.attrs = {
            'description': 'Normalized difference spectral index between bands at ' + str(prod.bvis) + ' and ' + str(
                prod.bswir) + ' nm', 'units': '-'}

        if allpixels:
            pass  # masked_raster = prod.raster.bands
        else:
            logging.info('apply water masking')
            mask = (ndwi_swir > prod.vis_swir_index_threshold) & (swir2 < prod.sunglint_threshold)  # (ndwi > -0.0) &
            masked = prod.raster.bands.where(mask)
            prod.raster['bands'] = masked
            prod.raster['sza'] = prod.raster['sza'].where(mask)

            # stop process if no valid (water) pixel available
            if np.isnan(prod.raster['sza'].values).all():
                logging.info('no water pixels, stop process')
                return

        ######################################
        # get DEM information
        # and adjust pressure raster
        ######################################
        if dem_file:
            logging.info('compute surface pressure from dem')
            dem = xr.open_dataset(dem_file).squeeze().interp(y=prod.raster.y,
                                                             x=prod.raster.x,
                                                             method='nearest')
            dem = dem.rename_vars({'band_data': 'dem'})
            dem.dem.attrs['long_name'] = 'digital elevation model'
            dem.dem.attrs['units'] = 'm'
            dem.dem.attrs['source'] = dem_file
            presure_msl = cams.raster.msl.interp(y=prod.raster.y, x=prod.raster.x)
            _pressure = (presure_msl * (1. - 0.0065 * dem.dem / 288.15) ** 5.255).values
        else:
            dem = None
            _pressure = cams.raster.sp.interp(x=prod.raster.x, y=prod.raster.y).values

        ######################################
        # load algorithm kernel
        ######################################
        kernel = Kernel(prod,
                        aero_lut,
                        trans_lut,
                        cams)

        if provide_kernel:
            return kernel

        ######################################
        # Aerosol type selection
        ######################################
        # select appropriate opac aerosol model from CAMS aod
        # remove URBAN for the moment
        opac_models = aero_lut.model.values
        models = aero_lut.drop_sel(model='URBA_rh70').model.values
        # get mean aot and aot550 from CAMS
        cams_aot_mean = cams.cams_aod.mean(['x', 'y'])
        cams_aot_ref = cams.cams_aod.interp(wl=550, method='quadratic')
        cams_aot_ref_mean = cams_aot_ref.mean(['x', 'y'])

        # get the model that has the closest aot spectral shape
        if opac_model is None:
            lut_aod = aero_lut.aot.sel(model=models, aot_ref=1).interp(wl=cams.cams_aod.wl)
            idx = np.abs((cams_aot_mean / cams_aot_ref_mean) - lut_aod).sum('wl').argmin()
            opac_model = aero_lut.sel(model=models).model.values[idx]

        logging.info('selected aerosol model: ' + opac_model)

        weights = (opac_models == opac_model).astype(int)

        ######################################
        # LUT preparation
        ######################################
        logging.info('lut interpolation')
        kernel.lut_preparation(weights=weights)
        kernel.set_gas_transmittance()

        ######################################
        # surface rugosity (wind) estimation
        ######################################
        # TODO speed up process (e.g. proper loop)
        logging.info('surface rugosity (wind) estimation')
        raster = kernel.get_coarse_masked_raster(xcoarsen=megapix_size,
                                                 ycoarsen=megapix_size)
        _Nwl, _height, _width = raster.bands.shape
        res = []
        for iy in range(0, _height, chunk):
            yc = min(_height, iy + chunk)

            for ix in range(0, _width, chunk):
                xc = min(_width, ix + chunk)
                res.append(kernel.rugosity_est_chunk(raster[dict(x=slice(ix, xc),
                                                                 y=slice(iy, yc))]))

        wind_img = xr.merge(res)
        # TODO add uncertainty weight for convolution
        kernel.wind_img = kernel.smoothing(wind_img,
                                           varname='wind',
                                           mask=np.ones((2, 2)),
                                           windows=np.array([1, 1]))


        ######################################
        # aerosol load estimation
        ######################################
        logging.info('aerosol optical thickness estimation')

        chunk = aot_chunk
        # kernel.aot_ref_cams_max=0.436
        raster = kernel.get_coarse_masked_raster(xcoarsen=aot_megapix_size, ycoarsen=aot_megapix_size)
        _Nwl, _height, _width = raster.bands.shape
        res = []
        for iy in range(0, _height, chunk):
            yc = min(_height, iy + chunk)

            for ix in range(0, _width, chunk):
                xc = min(_width, ix + chunk)
                raster_ = raster[dict(x=slice(ix, xc), y=slice(iy, yc))]
                kernel.aerosol_swir_chunk(raster_)
                #print('1', kernel.aot_ref_max)
                kernel.aerosol_swir_chunk(raster_, aot_refs=np.linspace(0, kernel.aot_ref_max, 21))
                #print('2', kernel.aot_ref_max)
                res.append(kernel.aerosol_visible_chunk(raster_))

        kernel.aot_ref_img = xr.merge(res)

        logging.info('smoothing aerosol optical thickness raster')
        # TODO add uncertainty weight for convolution
        aot_ref_raster = kernel.smoothing(kernel.aot_ref_img,
                                          varname='aot_ref',
                                          mask=np.ones((2*chunk,2*chunk)),
                                          windows=np.array([15, 15]))
        kernel.aot_ref_raster = aot_ref_raster

        ######################################
        # refined LUT preparation
        ######################################
        logging.info('refined LUT preparation')
        kernel.lut_preparation(weights=weights,
                               aot_refs=[0,*np.linspace(aot_ref_raster.aot_ref.min(),
                                                        aot_ref_raster.aot_ref.max(), 125)])

        ######################################
        # final grs process
        ######################################
        logging.info('final grs process')
        kernel.final_process(kernel.prod.raster,
                             aot_ref_raster.aot_ref)

        # filter wind retrievals
        kernel.wind_img['wind'] = kernel.wind_img.wind.where(kernel.wind_img.wind < 12)
        l2_prod = xr.merge([kernel.xres,
                  kernel.wind_img.interp(x=kernel.xres.x, y=kernel.xres.y)])

        ##############################################
        # Update flags and create mask from recipe
        ##############################################
        logging.info('Update flags and create mask from recipe')
        # flags for negative blue/green Rrs
        bitmask = 18
        prod.raster['flags'] = prod.raster.flags + (((l2_prod.Rrs.sel(wl=490, method='nearest') < -0.0005) |
                                                     (l2_prod.Rrs.sel(wl=565, method='nearest') < -0.0005)) << bitmask)
        # add name and description
        prod.raster.flags.attrs['flag_descriptions'][bitmask] = 'negative Rrs for blue or green bands'
        prod.raster.flags.attrs['flag_names'][bitmask] = 'neg_rrs'

        # mask from recipe
        mask = masking_.create_mask(prod.raster.flags,
                                    tomask=self.flags_tomask,
                                    tokeep=self.flags_tokeep,
                                    mask_name="mask")
        l2_prod = xr.merge([l2_prod, mask])

        l2_prod['central_wavelength'] = ('wl', prod.raster.wl_true.values)
        l2_prod = l2_prod.set_coords('central_wavelength')


        ######################################
        # Write final product
        ######################################
        logging.info('construct final product')
        # self.l2_prod = l2_prod
        self.l2a = L2aProduct(prod, l2_prod, cams, gas_trans, dem)
        del prod, l2_prod, cams, gas_trans, dem
        self.successful = True
        return kernel

    def write_output(self):
        logging.info('export final product into netcdf')
        self.l2a.export_to_netcdf(self.ofile,
                                  snap_compliant=self.snap_compliant)
