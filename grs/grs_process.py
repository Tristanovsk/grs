'''
Main program
'''

import numpy as np
import xarray as xr
import logging

from s2driver import driver_S2_SAFE as S2
from . import product, acutils, cams_product, l2a_product
from .fortran.grs import main_algo as grs_solver


class process:
    ''' '''

    def __init__(self):
        self.bandIds = range(13)
        self.type = np.float32
        self.chunk = 512

    def execute(self, file, ofile,
                cams_file='./cams.nc',
                resolution=20,
                allpixels=False,
                output='Rrs',
                logfile="log.txt",
                log_level="INFO",
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
        l1c = S2.s2image(file, band_idx=self.bandIds, resolution=resolution)
        l1c.load_product()
        logging.info('pass raw image as grs product object')
        prod = product(l1c.prod)

        ##################################
        # Fetch optional mask products
        ##################################
        # TODO (if necessary) activate MAJA reader and flag retrieval
        # prod.get_flags()

        ##################################
        # GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
        ##################################
        logging.info('get CAMS auxilliary data')
        cams = cams_product(prod, cams_file=cams_file)

        ##################################
        ## ADD ELEVATION AND PRESSURE BAND
        ##################################
        # TODO activate DEM loading to improve pressure computation
        # prod.get_elevation()

        #####################################
        # LOAD LUT FOR ATMOSPHERIC CORRECTION
        #####################################
        logging.info('loading look-up tables')
        lutf = acutils.lut(prod.band_names)
        lutc = acutils.lut(prod.band_names)
        lutf.load_lut(prod.lutfine, prod.sensordata.indband)
        lutc.load_lut(prod.lutcoarse, prod.sensordata.indband)

        ####################################
        #    absorbing gases correction
        ####################################
        logging.info('compute gaseous transmittance from cams data')
        gas_trans = acutils.gaseous_transmittance(prod, cams)
        Tg_raster = gas_trans.get_gaseous_transmittance()

        logging.info('correct for gaseous absorption')
        for wl in prod.wl.values:
            if (wl>1300) and (wl<1400):
                #do not correct outside atmospheric windows (e.g., cirrus bands)
                continue
            prod.raster['bands'].loc[wl] = prod.raster.bands.sel(wl=wl) / Tg_raster.sel(wl=wl).interp(x=prod.raster.x,
                                                                                                      y=prod.raster.y)
        prod.raster.bands.attrs['gas_absorption_correction'] = True

        ######################################
        # Water mask
        ######################################
        logging.info('compute spectral index (e.g., NDWI)')

        green = prod.raster.bands.sel(wl=prod.b565)
        nir = prod.raster.bands.sel(wl=prod.b865)
        swir = prod.raster.bands.sel(wl=prod.b1600)
        b2200 = prod.raster.bands.sel(wl=prod.b2200)

        ndwi = (green - nir) / (green + nir)
        ndwi_swir = (green - swir) / (green + swir)

        prod.raster['ndwi'] = ndwi
        prod.raster.ndwi.attrs = {
            'description': 'Normalized difference spectral index between bands at ' + str(prod.b565) + ' and ' + str(
                prod.b865) + ' nm', 'units': '-'}
        prod.raster['ndwi_swir'] = ndwi_swir
        prod.raster.ndwi_swir.attrs = {
            'description': 'Normalized difference spectral index between bands at ' + str(prod.b565) + ' and ' + str(
                prod.b1600) + ' nm', 'units': '-'}

        if allpixels:
            masked_raster = prod.raster.bands
        else:
            logging.info('apply water masking')
            prod.raster['bands'] = prod.raster.bands.where(ndwi > prod.ndwi_threshold). \
                where(b2200 < prod.sunglint_threshold). \
                where(ndwi_swir > prod.green_swir_index_threshold)
        self.raster = prod.raster

        ######################################
        # Rounding angle variables to save time in LUT interpolation
        ######################################
        logging.info('round angles for speed up lut interpolation')
        # def rounding(xarr, resol=1):
        #     vals = np.unique(xarr.round(resol))
        #     return vals[~np.isnan(vals)]
        #
        # sza_ = rounding(prod.raster.sza, 1)
        # azi_ = rounding((180 - prod.raster.raa) % 360, 0)
        # vza_ = rounding(prod.raster.vza, 1)
        wl_process = prod.wl_process
        Nband = len(wl_process)
        vza = prod.raster.vza.sel(wl=wl_process)
        sza = prod.raster.sza
        raa = prod.raster.raa.sel(wl=wl_process)

        sza_ = np.linspace(sza.min(), sza.max(), 10)
        vza_ = np.linspace(vza.min(), vza.max(), 20)
        raa_ = np.linspace(raa.min(), raa.max(), 60)

        ######################################
        # Set final parameters for grs processing
        ######################################
        logging.info('set final parameters for grs processing')

        eps_sunglint = prod.sensordata.rg
        rot = prod.sensordata.rot
        rrs = prod.rrs
        width = prod.width
        height = prod.height
        chunk = self.chunk

        logging.info('slice raster for desired wavelengths')
        raster = prod.raster['bands'].sel(wl=wl_process)

        solar_irr = prod.solar_irradiance.sel(wl=wl_process).values

        logging.info('get/set aerosol parameters')
        aotlut = np.array(lutf.aot, dtype=prod.type)
        fine_refl = lutf.refl.interp(vza=vza_).interp(azi=raa_).interp(sza=sza_)
        coarse_refl = lutc.refl.interp(vza=vza_).interp(azi=raa_).interp(sza=sza_)
        lut_shape = fine_refl.shape
        fine_Cext = lutf.Cext
        coarse_Cext = lutc.Cext
        aot_tot_cams_res = cams.cams_aod.interp(wavelength=wl_process)
        aot_sca_cams_res = aot_tot_cams_res * cams.cams_ssa.interp(wavelength=wl_process)

        # TODO implement pre-masking, now set to zero
        maskpixels = np.zeros((prod.height, prod.width))

        logging.info('get pressure full raster')

        ######################################
        # Run grs processing
        ######################################
        rcorr = np.full((Nband, width, height), np.nan, dtype=self.type)  # , order='F').T
        rcorrg = np.full((Nband, width, height), np.nan, dtype=self.type)  # , order='F').T
        aot550pix = np.full((width, height), np.nan, dtype=self.type)
        brdfpix = np.full((width, height), np.nan, dtype=self.type)

        logging.info('run grs process')
        for iy in range(0, width, chunk):
            yc = iy + chunk
            if yc > width:
                yc = width
            for ix in range(0, height, chunk):
                xc = ix + chunk
                if xc > height:
                    xc = height

                _sza = sza[ix:xc, iy:yc]
                nx, ny = _sza.shape
                if (nx == 0) or (ny == 0):
                    continue
                _raa = raa[:, ix:xc, iy:yc]
                _vza = vza[:, ix:xc, iy:yc]
                _maskpixels = maskpixels[ix:xc, iy:yc]
                _band_rad = raster[:, ix:xc, iy:yc]

                # prepare aerosol parameters
                aot_tot = aot_tot_cams_res.interp(x=_band_rad.x, y=_band_rad.y)
                aot_sca = aot_sca_cams_res.interp(x=_band_rad.x, y=_band_rad.y)
                aot550guess = cams.raster.aod550.interp(x=_band_rad.x, y=_band_rad.y)
                fcoef = np.full((nx, ny), 0.65)

                pressure_corr = cams.raster.sp.interp(x=_band_rad.x, y=_band_rad.y) * 1e-2 / prod.pressure_ref

                p = grs_solver.grs.main_algo(nx, ny, *lut_shape,
                                             aotlut, sza_, raa_, vza_,
                                             fine_refl, coarse_refl, fine_Cext, coarse_Cext,
                                             _vza, _sza, _raa, _band_rad.values, _maskpixels,
                                             wl_process, pressure_corr, eps_sunglint, solar_irr, rot,
                                             aot_tot, aot_sca, aot550guess, fcoef, rrs)

                rcorr[:, ix:xc, iy:yc] = p[0]
                rcorrg[:, ix:xc, iy:yc] = p[1]
                aot550pix[ix:xc, iy:yc] = p[2]
                brdfpix[ix:xc, iy:yc] = p[3]

        ######################################
        # Write final product
        ######################################
        logging.info('construct final product')

        Rrs = xr.DataArray(rcorr, coords=raster.coords, name='Rrs')
        Rrs_g = xr.DataArray(rcorrg, coords=raster.coords, name='Rrs_g')
        aot550 = xr.DataArray(aot550pix, coords={'y': raster.y, 'x': raster.x}, name='aot550')
        brdfg = xr.DataArray(brdfpix, coords={'y': raster.y, 'x': raster.x}, name='BRDFg')
        l2_prod = xr.merge([aot550, brdfg, Rrs, Rrs_g ])
        self.l2a = l2a_product(prod, l2_prod, cams, gas_trans)

        logging.info('export final product into netcdf')
        self.l2a.to_netcdf(ofile)
