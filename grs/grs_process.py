'''
Main program
'''

import numpy as np
import xarray as xr
import rioxarray as rio
import logging
import gc

from multiprocessing import Pool  # Process pool
from multiprocessing import sharedctypes
import itertools

from s2driver import driver_S2_SAFE as S2
# from .drivers import driver_S2_SAFE as S2
from . import product, acutils, auxdata, cams_product, l2a_product


# from .fortran.grs import main_algo as grs_solver


class process:
    ''' '''

    def __init__(self):
        self.bandIds = range(13)
        self.lut_file = '/data/vrtc/xlut/toa_lut_opac_wind_up.nc'
        self.Nproc = 38

    def execute(self, file, ofile,
                cams_file='./cams.nc',
                surfwater_file=None,
                resolution=20,
                allpixels=False,
                snap_compliant=False
                ):

        '''
        Main program calling all GRS steps

        :param file: Input file to be processed
        :param ofile: Absolute path of the output file
        :param cams_file: Absolute path for root directory of CAMS data

        :param sensor: Set the sensor type: S2A, S2B, LANDSAT_5, LANDSAT_7, LANDSAT_8
                    (by default sensor type is retrieved from input file name)
        :param resolution: pixel resolution in meters (integer), choice between: 10, 20, 60 m
        :param allpixels: force to process all pixels even they are flagged as "Vegetation" or "Non-water"
        :param output: set the unit of the retrievals:
                 * 'Lwn', normalized water-leaving radiance (in  :math:`mW cm^{-2} sr^{-1} \mu m^{-1})`
                 * 'Rrs', remote sensing reflectance (in  :math:`sr^{-1}`)
                 {default: 'Rrs'}
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
        self.prod=prod
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
        cams = cams_product(prod.raster, cams_file=cams_file)

        # get aerosol optical thickness at 550nm
        cams_aot_ref = cams.cams_aod.interp(wl=550, method='quadratic').compute()

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
        aero_lut = xr.open_dataset(self.lut_file)
        aero_lut['wl'] = aero_lut['wl'] * 1000
        aero_lut['aot'] = aero_lut.aot.isel(wind=0).squeeze()
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

        #select appropriate opac aerosol model from CAMS aod
        #remove URBAN for the moment
        models = aero_lut.drop_sel(model='URBA_rh70').model.values
        # get mean aot and aot550 from CAMS
        cams_aot_mean = cams.cams_aod.mean(['x', 'y'])
        cams_aot_ref = cams.cams_aod.interp(wl=550, method='quadratic')
        cams_aot_ref_mean = cams_aot_ref.mean(['x', 'y'])

        # get the model that has the closest aot spactral shape
        lut_aod = aero_lut.aot.sel(model=models, aot_ref=1).interp(wl=cams.cams_aod.wl)
        idx = np.abs((cams_aot_mean / cams_aot_ref_mean) - lut_aod).sum('wl').argmin()
        opac_model = aero_lut.sel(model=models).model.values[idx]

        # slice LUT
        aero_lut_ = aero_lut.sel(wind=_wind, method='nearest').sel(model=opac_model)

        # get AOT550 raster (TODO replace with optimal estimation)
        aot_ref_raster = cams_aot_ref
        ang_resol = {'sza': 1, 'vza': 1, 'raa': 0}
        for param in ['sza', 'vza', 'raa']:
            prod.raster[param + '_trunc'] = prod.raster[param].round(ang_resol[param])

        sza_ = np.unique(prod.raster.sza_trunc)
        vza_ = np.unique(prod.raster.vza_trunc)
        azi_ = (180 - np.unique(prod.raster.raa_trunc)) % 360
        sza_ = sza_[~np.isnan(sza_)]
        vza_ = vza_[~np.isnan(vza_)]
        azi_ = azi_[~np.isnan(azi_)]

        sza_lut_step = 2
        vza_lut_step = 2
        sza_slice = slice(np.min(sza_) - sza_lut_step, np.max(sza_) + sza_lut_step)
        vza_slice = slice(np.min(vza_) - vza_lut_step, np.max(vza_) + vza_lut_step)
        tweak = 4
        aot_ref_ = np.unique((aot_ref_raster / tweak).round(3)) * tweak
        aot_ref_min = aot_ref_raster.min()
        aot_ref_max = aot_ref_raster.max()
        aot_lut = aero_lut_.aot.interp(wl=wl_true, method='quadratic')
        aot_lut = aot_lut.interp(aot_ref=np.linspace(aot_ref_min, aot_ref_max, 1000))  # .plot(hue='wl')

        Rdiff_lut = aero_lut_.I.sel(sza=sza_slice, vza=vza_slice).interp(wl=wl_true, method='quadratic')
        Rdiff_lut = Rdiff_lut.interp(sza=sza_).interp(aot_ref=aot_ref_, method='quadratic').interp(vza=vza_).interp(
            azi=azi_)

        szas = Rdiff_lut.sza.values
        vzas = Rdiff_lut.vza.values
        azis = Rdiff_lut.azi.values
        aot_refs = Rdiff_lut.aot_ref.values

        ######################################
        # Set final parameters for grs processing
        ######################################

        width = prod.width
        height = prod.height
        Nwl = len(prod.raster.wl_to_process)

        pressure_ref = 101500.

        _sunglint_eps = sunglint_eps.values

        # prepare aerosol parameters
        aot_ref_raster = cams_aot_ref.interp(x=prod.raster.x, y=prod.raster.y).drop('wl')
        aot_ref_raster.name = 'aot550'
        # _aot = aot_lut.interp(aot_ref=_aot_ref)
        _pressure = cams.raster.sp.interp(x=prod.raster.x, y=prod.raster.y).values
        _rot = rot.values

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
            arr_tmp = np.full((Nwl, Ny, Nx), np.nan, dtype=self.prod._type)

            # subsetting
            _sza = prod.raster.sza[iy:yc, ix:xc]  # .values
            _raa = prod.raster.raa[:, iy:yc, ix:xc]
            _azi = (180. - _raa) % 360
            _vza = prod.raster.vza[:, iy:yc, ix:xc]
            _air_mass_ = acutils.misc.air_mass(_sza,
                                               _vza).values  # air_mass[:, iy:yc,ix:xc] #air_mass(_sza,_vza).values #_p_slope = prod.raster.p_slope[:, iy:yc,ix:xc]
            _p_slope_ = prod.p_slope(_sza, _vza, _raa, sigma2=_sigma2).values  # _p_slope[:, iy:yc,ix:xc]
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
            _aot = acutils._interp_aotlut(aot_lut.aot_ref.values, _aot_ref, Nwl, Ny, Nx, aot_lut.values)

            #  correction for diffuse light
            Rcorr = _band_rad.values - _Rdiff

            # direct transmittance up/down
            Tdir = acutils.misc.transmittance_dir(_aot, _air_mass_, _rot_raster)

            Rf = np.full((len(self.prod.iwl_swir), Ny, Nx), np.nan, dtype=self.prod._type)

            for iwl in self.prod.iwl_swir:
                Rf[iwl] = Rcorr[iwl] / (Tdir[iwl] * _sunglint_eps[iwl] * _p_slope_[iwl])
            Rf[Rf < 0] = 0.
            Rf = np.min(Rf, axis=0)
            Rf_tmp[iy:yc, ix:xc] = Rf
            Rf = acutils._multiplicate(_sunglint_eps, Rf, arr_tmp)
            Rf = Tdir * _p_slope_ * Rf

            Rrs_tmp[:, iy:yc, ix:xc] = ((Rcorr - Rf) / np.pi)
            return

        window_idxs = [(i, j) for i, j in
                       itertools.product(range(0, width, self.prod.chunk),
                                         range(0, height, self.prod.chunk))]

        global pool
        pool = Pool(self.Nproc)
        res = pool.map(chunk_process, window_idxs)
        pool.terminate()
        pool = None
        print('success')

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
