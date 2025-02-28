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

from . import Product, acutils, AuxData, CamsProduct, L2aProduct, Masking, Rasterization

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
                scale_aot=1,
                opac_model=None,
                allpixels=False,
                snap_compliant=False
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
                prod = Product(xr.open_dataset(l1c_prod), engine=NETCDF_ENGINE)
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

        # Cox-Munk isotropic mean square slope (sigma2)
        wind = np.sqrt(cams.raster['v10'] ** 2 + cams.raster['u10'] ** 2) * 2
        sigma2 = (wind + 0.586) / 195.3

        # get mean values to set LUT
        _sigma2 = sigma2.median().values
        _wind = wind.mean().values

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
        Ttot_Ed = xr.open_dataset(self.trans_lut_file, engine=NETCDF_ENGINE)
        Ttot_Ed['wl'] = Ttot_Ed['wl'] * 1000

        aero_lut = xr.open_dataset(self.lut_file, engine=NETCDF_ENGINE)
        aero_lut['wl'] = aero_lut['wl'] * 1000
        aero_lut['aot'] = aero_lut.aot.isel(wind=0).squeeze()
        aero_lut['aaot'] = aero_lut['aot'] * (1 - aero_lut.ssa.isel(wind=0).squeeze())

        # remove URBAN aerosol model for this example
        models = aero_lut.drop_sel(model='URBA_rh70').model.values

        _auxdata = AuxData(wl=wl_true)  # wl=masked.wl)
        sunglint_eps = _auxdata.sunglint_eps  # ['mean'].interp(wl=wl_true)
        rot = _auxdata.rot

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
        # LUT preparation
        ######################################
        logging.info('lut interpolation')

        # select appropriate opac aerosol model from CAMS aod
        # remove URBAN for the moment
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
        # slice LUT
        aero_lut_ = aero_lut.sel(wind=_wind, method='nearest').sel(model=opac_model)

        # get AOT550 raster (TODO replace with optimal estimation)
        logging.info('scaling aot by: ' + str(scale_aot))
        aot_ref_raster = cams_aot_ref * scale_aot
        aot_ref_raster = aot_ref_raster.astype(np.float32)

        # get unique values for angles and further lut interpolation
        ang_resol = {'sza': 0.1, 'vza': 0.1, 'raa_round': 0}
        szamin, szamax = float(prod.raster['sza'].min()), float(prod.raster['sza'].max())
        vzamin, vzamax = float(prod.raster.isel(wl=0)['vza'].min()), float(prod.raster.isel(wl=0)['vza'].max())

        # check for out-of-range
        def check_out_of_range(vmin, vmax, ceiling=88):
            vmin = np.max([0, vmin])
            vmax = np.min([ceiling, vmax])
            return vmin, vmax

        szamin, szamax = check_out_of_range(szamin, szamax)
        vzamin, vzamax = check_out_of_range(vzamin, vzamax, ceiling=25)

        sza_ = np.arange(szamin, szamax + ang_resol['sza'], ang_resol['sza'])
        vza_ = np.arange(vzamin, vzamax + ang_resol['vza'], ang_resol['vza'])

        azi_ = (180 - np.unique(prod.raster.isel(wl=0)['raa'].round(ang_resol['raa_round']))) % 360
        azi_ = azi_[~np.isnan(azi_)]

        sza_lut_step = 2
        vza_lut_step = 2

        sza_slice = slice(np.min(sza_) - sza_lut_step, np.max(sza_) + sza_lut_step)
        vza_slice = slice(np.min(vza_) - vza_lut_step, np.max(vza_) + vza_lut_step)

        tweak = 4
        aot_ref_ = np.unique((aot_ref_raster / tweak).round(3)) * tweak
        aot_ref_min = 0.  # aot_ref_raster.min()
        aot_ref_max = aot_ref_raster.max()
        aot_lut = aero_lut_.aot.interp(wl=wl_true, method='quadratic')
        aot_lut = aot_lut.interp(aot_ref=np.linspace(aot_ref_min, aot_ref_max.values, 1000))  # .plot(hue='wl')
        aaot_lut = aero_lut_.aaot.interp(wl=wl_true, method='quadratic')
        aaot_lut = aaot_lut.interp(aot_ref=np.linspace(aot_ref_min, aot_ref_max.values, 1000))

        Rdiff_lut = aero_lut_.I.sel(sza=sza_slice,
                                    vza=vza_slice
                                    ).interp(wl=wl_true,
                                             method='quadratic'
                                             ).interp(azi=azi_)
        Rdiff_lut = Rdiff_lut.interp(sza=sza_, vza=vza_)
        Rray = Rdiff_lut.sel(aot_ref=0)
        Rdiff_lut = Rdiff_lut.interp(aot_ref=[0, 0.02, 0.05, 0.07, *aot_ref_],
                                     method='quadratic') #.sortby("aot_ref")

        szas = Rdiff_lut.sza.values
        vzas = Rdiff_lut.vza.values
        azis = Rdiff_lut.azi.values
        aot_refs = Rdiff_lut.aot_ref.values

        Ttot_Ed_ = Ttot_Ed.sel(model=opac_model).sel(wind=_wind, method='nearest').interp(sza=szas).interp(
            aot_ref=aot_ref_, method='quadratic').interp(wl=wl_true, method='cubic').Ttot_Ed
        Ttot_Lu_ = Ttot_Ed.sel(model=opac_model).sel(wind=_wind, method='nearest').interp(sza=vzas).interp(
            aot_ref=aot_ref_, method='quadratic').interp(wl=wl_true, method='cubic').Ttot_Ed ** 1.05

        ######################################
        # Set final parameters for grs processing
        ######################################
        logging.info('set final parameters')
        width = prod.width
        height = prod.height
        Nwl = len(prod.raster.wl_to_process)

        pressure_ref = self.pressure_ref

        _sunglint_eps = sunglint_eps.values

        # prepare aerosol parameters
        aot_ref_raster = aot_ref_raster.interp(x=prod.raster.x,
                                               y=prod.raster.y).drop('wl').astype(np.float32)
        aot_ref_raster.name = 'aot550'
        _rot = rot.values

        # _aot = aot_lut.interp(aot_ref=_aot_ref)
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

        # -------------------------------------------------------------
        # SET GASEOUS TRANSMITTANCE FOR LOW ALTITUDE GASES
        # -------------------------------------------------------------
        gases = ['h2o', 'ch4']
        gas_trans = acutils.GaseousTransmittance(prod, cams)
        # set total transmittance values
        gas_trans.coef_abs_scat['h2o'] = 0.5
        gas_trans.coef_abs_scat['ch4'] = 0.65
        Tg_diff_raster = gas_trans.get_gaseous_transmittance(gases=gases, background=False).transpose("wl", "y", "x")
        Tg_diff_raster = Tg_diff_raster.interp(x=prod.raster.x, y=prod.raster.y)

        # set total transmittance values
        for gas in gases:
            gas_trans.coef_abs_scat[gas] = 1
        Tg_raster = gas_trans.get_gaseous_transmittance(gases=gases, background=False).transpose("wl", "y", "x")
        Tg_raster = Tg_raster.interp(x=prod.raster.x, y=prod.raster.y)

        ######################################
        # Run grs processing
        ######################################
        logging.info('run grs process')
        global chunk_process
        Rrs_result = np.ctypeslib.as_ctypes(np.full((Nwl, height, width), np.nan, dtype=prod._type))
        Rf_result = np.ctypeslib.as_ctypes(np.full((height, width), np.nan, dtype=prod._type))
        shared_Rrs = sharedctypes.RawArray(Rrs_result._type_, Rrs_result)
        shared_Rf = sharedctypes.RawArray(Rf_result._type_, Rf_result)

        def chunk_process(args):
            iy, ix = args
            yc = min(height, iy + prod.chunk)
            xc = min(width, ix + prod.chunk)
            Rrs_tmp = np.ctypeslib.as_array(shared_Rrs)
            Rf_tmp = np.ctypeslib.as_array(shared_Rf)

            _band_rad = prod.raster.bands[:, iy:yc, ix:xc]

            Nwl, Ny, Nx = _band_rad.shape
            if Ny == 0 or Nx == 0:
                return
            arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=prod._type)

            # subsetting
            _sza = prod.raster.sza[iy:yc, ix:xc]  # .values
            if monoview:
                _raa = prod.raster.raa[iy:yc, ix:xc]
                _vza = prod.raster.vza[iy:yc, ix:xc]
                _vza_mean = _vza.values
            else:
                _raa = prod.raster.raa[:, iy:yc, ix:xc]
                _vza = prod.raster.vza[:, iy:yc, ix:xc]
                _vza_mean = np.mean(_vza, axis=0).values

            _azi = (180. - _raa) % 360
            _air_mass_ = acutils.Misc.air_mass(_sza, _vza).values
            _p_slope_ = prod.p_slope(_sza, _vza, _raa, sigma2=_sigma2, monoview=monoview).values
            _aot_ref = aot_ref_raster.values[iy:yc, ix:xc]
            _pressure_ = _pressure[iy:yc, ix:xc] / pressure_ref
            _Tg_abs = Tg_raster[:, iy:yc, ix:xc].values
            _Tg_abs_diff = Tg_diff_raster[:, iy:yc, ix:xc].values

            # construct wl,y,x raster for Rayleigh optical thickness
            _rot_raster = _R_._multiplicate(_rot, _pressure_, arr_tmp)

            # get LUT values
            _Rdiff = _R_.interp_Rlut(szas, _sza.values,
                                     vzas, _vza.values,
                                     azis, _azi.values,
                                     aot_refs, _aot_ref,
                                     Nwl, Ny, Nx, Rdiff_lut.values)

            _Rray = _R_.interp_Rlut_rayleigh(szas, _sza.values,
                                             vzas, _vza.values,
                                             azis, _azi.values,
                                             Nwl, Ny, Nx, Rray.values)

            # _Rdiff = _Rdiff + (_pressure_ - 1) * _Rray
            _Rdiff = _Rdiff * _Tg_abs_diff * _pressure_

            _aot = _R_._interp_aotlut(aot_lut.aot_ref.values, _aot_ref, Nwl, Ny, Nx, aot_lut.values)

            #  correction for diffuse light
            Rcorr = _band_rad.values - _Rdiff

            # direct transmittance up/down
            Tdir = acutils.Misc.transmittance_dir(_aot, _air_mass_, _rot_raster)

            # vTotal transmittance (for Ed and Lu)
            Tdown = _R_._interp_Tlut(szas, _sza.values, Ttot_Ed_.aot_ref.values, _aot_ref, Nwl, Ny, Nx,
                                     Ttot_Ed_.values)
            Tup = _R_._interp_Tlut(vzas, _vza_mean, Ttot_Ed_.aot_ref.values, _aot_ref, Nwl, Ny, Nx, Ttot_Lu_.values)
            Ttot_du = Tdown * Tup * _Tg_abs

            Rf = np.full((len(prod.iwl_swir), Ny, Nx), np.nan, dtype=prod._type)

            for iwl in prod.iwl_swir:

                if monoview:
                    Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _Tg_abs[iwl] * _sunglint_eps[iwl] * _p_slope_)
                else:
                    Rf[iwl] = (_sunglint_eps[-1] * _p_slope_[-1] * Rcorr[iwl] /
                                (Tdir[iwl] * _Tg_abs[iwl] * _sunglint_eps[iwl] * _p_slope_[iwl]))

            Rf[Rf < 0] = 0.
            Rf = np.min(Rf, axis=0)
            Rf_tmp[iy:yc, ix:xc] = Rf

            Rf = _R_._multiplicate(_sunglint_eps, Rf, arr_tmp)
            Rf = _Tg_abs * Tdir * Rf * _p_slope_ / (_sunglint_eps[-1] * _p_slope_[-1])

            # sunglint removal
            #Rrs_tmp_ =Rcorr / np.pi# Rrs_tmp[:, iy:yc, ix:xc]
            Rrs_tmp_ = ((Rcorr - Rf) / np.pi)

            # Convert from TOA to BOA for positive values
            Ttot_du[Rrs_tmp_ < 0] = 1.
            Rrs_tmp_ = Rrs_tmp_ / Ttot_du
            Rrs_tmp[:, iy:yc, ix:xc] = Rrs_tmp_
            return

        window_idxs = [(i, j) for i, j in
                       itertools.product(range(0, height, prod.chunk),
                                         range(0, width, prod.chunk))]

        global pool
        pool = Pool(self.Nproc)
        res = pool.map(chunk_process, window_idxs)
        pool.terminate()
        pool = None
        logging.info('success')

        ######################################
        # construct l2a object
        ######################################
        logging.info('construct final product')
        self.aot_ref_raster = aot_ref_raster
        l2_prod = xr.Dataset(dict(Rrs=(['wl', "y", "x"], np.ctypeslib.as_array(shared_Rrs)),
                                  BRDFg=(["y", "x"], np.ctypeslib.as_array(shared_Rf)),
                                  aot550=(["y", "x"], aot_ref_raster.values)),
                             coords=dict(wl=prod.raster.wl,
                                         x=prod.raster.x,
                                         y=prod.raster.y),
                             )

        l2_prod['central_wavelength'] = ('wl', prod.raster.wl_true.values)
        l2_prod = l2_prod.set_coords('central_wavelength')

        ##############################################
        # Update flags and create mask from recipe
        ##############################################
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

        ######################################
        # Write final product
        ######################################
        logging.info('construct final product')
        # self.l2_prod = l2_prod
        self.l2a = L2aProduct(prod, l2_prod, cams, gas_trans, dem)
        del prod, l2_prod, cams, gas_trans, dem
        self.successful = True
        return

    def write_output(self):
        logging.info('export final product into netcdf')
        self.l2a.export_to_netcdf(self.ofile,
                                  snap_compliant=self.snap_compliant)
