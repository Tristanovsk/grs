'''
Main program
'''

import os

import importlib_resources
import yaml

import numpy as np
import xarray as xr
import rioxarray as rio
import logging
import gc

from multiprocessing import Pool  # Process pool
from multiprocessing import sharedctypes
import itertools

from s2driver import driver_S2_SAFE as S2

from . import product, acutils, auxdata, cams_product, l2a_product

opj = os.path.join

configfile = importlib_resources.files(__package__) / 'config.yml'
with open(configfile, 'r') as file:
    config = yaml.safe_load(file)

GRSDATA = config['path']['grsdata']
TOALUT = config['path']['toa_lut']
TRANSLUT = config['path']['trans_lut']
CAMS_PATH = config['path']['trans_lut']
NCPU = config['processor']['ncpu']


class process:
    ''' '''

    def __init__(self):
        self.bandIds = range(13)
        self.lut_file = opj(GRSDATA, TOALUT)
        self.trans_lut_file = opj(GRSDATA, TRANSLUT)
        self.cams_dir = CAMS_PATH
        self.Nproc = NCPU

    def execute(self, file, ofile,
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
        :param file:  Input file to be processed
        :param ofile: Absolute path of the output file
        :param cams_file: Absolute path for root directory of CAMS data
        :param surfwater_file: Absolute path the surfwater filr (.tif)
        :param dem_file: Absolute path of the DEM geotiff file
        :param resolution: pixel resolution in meter
        :param scale_aot: scaling factor applied to CAMS aod550 raster
        :param opac_model: If set force OPAC aerosol model for LUT interpolation (taken from CAMS data otherwise)
        :param allpixels: if True process all pixels (no water pixel masking)
        :param snap_compliant: Output format compliant with SNAP software for practical analysis
        :return:
        '''



        ##################################
        # Get image data
        ##################################
        logging.info('Open raw image and compute angle parameters')
        global l1c
        l1c = S2.sentinel2_driver(file, band_idx=self.bandIds, resolution=resolution)
        l1c.load_product()
        logging.info('pass raw image as grs product object')
        prod = product(l1c.prod)
        self.prod = prod
        wl_true = prod.raster.wl_true

        # clear memory (TODO make it work!!)
        del l1c
        gc.collect()

        ##################################
        # Fetch optional mask products
        ##################################
        # TODO (if necessary) activate MAJA reader and flag retrieval
        # prod.get_flags()

        ##################################
        # GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
        ##################################
        logging.info('get CAMS auxilliary data')
        if cams_file:
            cams = cams_product(prod.raster, cams_file=cams_file)
        else:
            tile = prod.raster.attrs['tile']
            cams_dir = os.path.join(self.cams_dir, tile)
            cams = cams_product(prod.raster, dir=cams_dir, suffix='_' + tile)
        cams.load()

        # Cox-Munk isotropic mean square slope (sigma2)
        wind = np.sqrt(cams.raster['v10'] ** 2 + cams.raster['u10'] ** 2)
        sigma2 = (wind + 0.586) / 195.3

        # get mean values to set LUT
        _sigma2 = sigma2.mean().values
        _wind = wind.mean().values

        if surfwater_file:
            prod.raster.surfwater = rio.open_rasterio(surfwater_file).squeeze()
            prod.raster.surfwater = prod.raster.surfwater.where(prod.raster.surfwater != 255., 32).astype(np.int8)
            prod.raster.surfwater = prod.raster.surfwater.interp(x=prod.x, y=prod.y, method='nearest')
            prod.raster.surfwater.name = 'surfwater'
            prod.raster.surfwater.attrs = {
                'description': 'surfwater file not provided as input, all pixels flagged as water (e.g., surfwater=1)'}

        #####################################
        # SUBSET RASTER TO KEEP REQUESTED BANDS
        #####################################
        if prod.bcirrus:
            prod.cirrus = prod.raster.bands.sel(wl=prod.bcirrus, method='nearest')
        if prod.bwv:
            prod.wv = prod.raster.bands.sel(wl=prod.bwv, method='nearest')

        prod.raster = prod.raster.sel(wl=prod.wl_process, method='nearest')

        ##################################
        ## ADD ELEVATION AND PRESSURE BAND
        ##################################
        # TODO activate DEM loading to improve pressure computation
        # prod.get_elevation()

        #####################################
        # LOAD LUT FOR ATMOSPHERIC CORRECTION
        #####################################
        logging.info('loading look-up tables')
        Ttot_Ed = xr.open_dataset(self.trans_lut_file)
        Ttot_Ed['wl'] = Ttot_Ed['wl'] * 1000

        aero_lut = xr.open_dataset(self.lut_file)
        aero_lut['wl'] = aero_lut['wl'] * 1000
        aero_lut['aot'] = aero_lut.aot.isel(wind=0).squeeze()

        # remove URBAN aerosol model for this example
        models = aero_lut.drop_sel(model='URBA_rh70').model.values

        _auxdata = auxdata(wl=wl_true)  # wl=masked.wl)
        sunglint_eps = _auxdata.sunglint_eps  # ['mean'].interp(wl=wl_true)
        rot = _auxdata.rot

        ####################################
        #    absorbing gases correction
        ####################################
        logging.info('compute gaseous transmittance from cams data')
        gas_trans = acutils.gaseous_transmittance(prod, cams)
        Tg_raster = gas_trans.get_gaseous_transmittance()

        logging.info('correct for gaseous absorption')
        for wl in prod.raster.wl.values:
            prod.raster['bands'].loc[wl] = prod.raster.bands.sel(wl=wl) / Tg_raster.sel(wl=wl).interp(x=prod.raster.x,
                                                                                                      y=prod.raster.y)
        prod.raster.bands.attrs['gas_absorption_correction'] = True

        ######################################
        # Water mask
        ######################################
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

        self.raster = prod.raster

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

        logging.info('selected aerosol model: '+opac_model)
        # slice LUT
        aero_lut_ = aero_lut.sel(wind=_wind, method='nearest').sel(model=opac_model)

        # get AOT550 raster (TODO replace with optimal estimation)
        logging.info('scaling aot by: ' + str(scale_aot))
        aot_ref_raster = cams_aot_ref * scale_aot

        # get unique values for angles and further lut interpolation
        ang_resol = {'sza': 0.1, 'vza': 0.1, 'raa_round': 0}
        szamin, szamax = float(prod.raster['sza'].min()), float(prod.raster['sza'].max())
        vzamin, vzamax = float(prod.raster['vza'].isel(wl=0).min()), float(prod.raster['vza'].isel(wl=0).max())

        # check for out-of-range
        def check_out_of_range(vmin, vmax, ceiling=88):
            vmin = np.max([0, vmin])
            vmax = np.min([ceiling, vmax])
            return vmin, vmax

        szamin, szamax = check_out_of_range(szamin, szamax)
        vzamin, vzamax = check_out_of_range(vzamin, vzamax, ceiling=25)

        sza_ = np.arange(szamin, szamax + ang_resol['sza'], ang_resol['sza'])
        vza_ = np.arange(vzamin, vzamax + ang_resol['vza'], ang_resol['vza'])

        azi_ = (180 - np.unique(prod.raster['raa'].isel(wl=0).round(ang_resol['raa_round']))) % 360
        azi_ = azi_[~np.isnan(azi_)]

        sza_lut_step = 2
        vza_lut_step = 2

        sza_slice = slice(np.min(sza_) - sza_lut_step, np.max(sza_) + sza_lut_step)
        vza_slice = slice(np.min(vza_) - vza_lut_step, np.max(vza_) + vza_lut_step)

        tweak = 3
        aot_ref_ = np.unique((aot_ref_raster / tweak).round(3)) * tweak
        aot_ref_min = aot_ref_raster.min()
        aot_ref_max = aot_ref_raster.max()
        aot_lut = aero_lut_.aot.interp(wl=wl_true, method='quadratic')
        aot_lut = aot_lut.interp(aot_ref=np.linspace(aot_ref_min, aot_ref_max, 1000))  # .plot(hue='wl')

        Rdiff_lut = aero_lut_.I.sel(sza=sza_slice, vza=vza_slice).interp(wl=wl_true, method='quadratic').interp(
            azi=azi_)
        Rdiff_lut = Rdiff_lut.interp(sza=sza_, vza=vza_)
        Rray = Rdiff_lut.sel(aot_ref=0)
        Rdiff_lut = Rdiff_lut.interp(aot_ref=aot_ref_, method='quadratic')

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

        pressure_ref = 101500.

        _sunglint_eps = sunglint_eps.values

        # prepare aerosol parameters
        aot_ref_raster = aot_ref_raster.interp(x=prod.raster.x, y=prod.raster.y).drop('wl')
        aot_ref_raster.name = 'aot550'
        _rot = rot.values

        # _aot = aot_lut.interp(aot_ref=_aot_ref)
        if dem_file:
            logging.info('compute surface pressure from dem')
            dem = xr.open_dataset(dem_file).squeeze().interp(y=prod.raster.y, x=prod.raster.x, method='nearest')
            dem = dem.rename_vars({'band_data': 'elevation'})
            presure_msl = cams.raster.msl.interp(y=prod.raster.y, x=prod.raster.x)
            _pressure = (presure_msl * (1. - 0.0065 * dem.elevation / 288.15) ** 5.255).values
        else:
            _pressure = cams.raster.sp.interp(x=prod.raster.x, y=prod.raster.y).values

        ######################################
        # Run grs processing
        ######################################
        logging.info('run grs process')
        global chunk_process
        Rrs_result = np.ctypeslib.as_ctypes(np.full((Nwl, width, height), np.nan, dtype=self.prod._type))
        Rf_result = np.ctypeslib.as_ctypes(np.full((width, height), np.nan, dtype=self.prod._type))
        shared_Rrs = sharedctypes.RawArray(Rrs_result._type_, Rrs_result)
        shared_Rf = sharedctypes.RawArray(Rf_result._type_, Rf_result)

        def chunk_process(args):
            iy, ix = args
            yc = min(height, iy + self.prod.chunk)
            xc = min(width, ix + self.prod.chunk)
            Rrs_tmp = np.ctypeslib.as_array(shared_Rrs)
            Rf_tmp = np.ctypeslib.as_array(shared_Rf)

            _band_rad = prod.raster.bands[:, iy:yc, ix:xc]

            Nwl, Ny, Nx = _band_rad.shape
            if Ny == 0 or Nx == 0:
                return
            arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=self.prod._type)

            # subsetting
            _sza = prod.raster.sza[iy:yc, ix:xc]  # .values
            _raa = prod.raster.raa[:, iy:yc, ix:xc]
            _azi = (180. - _raa) % 360
            _vza = prod.raster.vza[:, iy:yc, ix:xc]
            _vza_mean = np.mean(_vza, axis=0).values
            _air_mass_ = acutils.misc.air_mass(_sza, _vza).values
            _p_slope_ = prod.p_slope(_sza, _vza, _raa, sigma2=_sigma2).values
            _aot_ref = aot_ref_raster.values[iy:yc, ix:xc]
            _pressure_ = _pressure[iy:yc, ix:xc] / pressure_ref

            # construct wl,y,x raster for Rayleigh optical thickness
            _rot_raster = acutils._multiplicate(_rot, _pressure_, arr_tmp)

            # get LUT values
            _Rdiff = acutils._interp_Rlut(szas, _sza.values,
                                          vzas, _vza.values,
                                          azis, _azi.values,
                                          aot_refs, _aot_ref,
                                          Nwl, Ny, Nx, Rdiff_lut.values)

            _Rray = acutils._interp_Rlut_rayleigh(szas, _sza.values,
                                                  vzas, _vza.values,
                                                  azis, _azi.values,
                                                  Nwl, Ny, Nx, Rray.values)

            _Rdiff = _Rdiff + (_pressure_ - 1) * _Rray

            _aot = acutils._interp_aotlut(aot_lut.aot_ref.values, _aot_ref, Nwl, Ny, Nx, aot_lut.values)

            #  correction for diffuse light
            Rcorr = _band_rad.values - _Rdiff

            # direct transmittance up/down
            Tdir = acutils.misc.transmittance_dir(_aot, _air_mass_, _rot_raster)

            # vTotal transmittance (for Ed and Lu)
            Tdown = acutils._interp_Tlut(szas, _sza.values, Ttot_Ed_.aot_ref.values, _aot_ref, Nwl, Ny, Nx,
                                         Ttot_Ed_.values)
            Tup = acutils._interp_Tlut(vzas, _vza_mean, Ttot_Ed_.aot_ref.values, _aot_ref, Nwl, Ny, Nx, Ttot_Lu_.values)
            Ttot_du = Tdown * Tup

            Rf = np.full((len(self.prod.iwl_swir), Ny, Nx), np.nan, dtype=self.prod._type)

            for iwl in self.prod.iwl_swir:
                Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _sunglint_eps[iwl] * _p_slope_[iwl])
            Rf[Rf < 0] = 0.
            Rf = np.min(Rf, axis=0)
            Rf_tmp[iy:yc, ix:xc] = Rf
            Rf = acutils._multiplicate(_sunglint_eps, Rf, arr_tmp)
            Rf = Tdir * _p_slope_ * Rf

            Rrs_tmp[:, iy:yc, ix:xc] = ((Rcorr - Rf) / np.pi) / Ttot_du
            return

        window_idxs = [(i, j) for i, j in
                       itertools.product(range(0, width, self.prod.chunk),
                                         range(0, height, self.prod.chunk))]

        global pool
        pool = Pool(self.Nproc)
        res = pool.map(chunk_process, window_idxs)
        pool.terminate()
        pool = None
        logging.info('success')

        ######################################
        # Write final product
        ######################################
        logging.info('construct final product')
        xres = xr.Dataset(dict(Rrs=(['wl', "y", "x"], np.ctypeslib.as_array(shared_Rrs)),
                               BRDFg=(["y", "x"], np.ctypeslib.as_array(shared_Rf))),
                          coords=dict(wl=prod.raster.wl,
                                      x=prod.raster.x,
                                      y=prod.raster.y),
                          )

        l2_prod = xr.merge([aot_ref_raster, xres])
        self.l2a = l2a_product(prod, l2_prod, cams, gas_trans)

        logging.info('export final product into netcdf')
        self.l2a.to_netcdf(ofile, snap_compliant=snap_compliant)

        return
