import os
import xarray as xr
import rioxarray as rio
import logging


class l2a_product():
    def __init__(self, prod, l2_prod, cams, gas_trans):
        self.prod = prod
        self.l2_prod = l2_prod
        self.cams = cams
        self.gas_trans = gas_trans
        self.wl_cirrus_band = 1375
        self.wl_o2_band = 945
        self.complevel = 5
        self.construct_l2a()

    def construct_l2a(self):
        logging.info('construct l2a')
        native_raster = self.prod.raster
        # first get attributes from native product
        self.l2_prod.attrs = native_raster.attrs

        keys = ['sunglint_threshold', 'ndwi_threshold', 'green_swir_index_threshold', 'hcld_threshold',
                'dirdata', 'abs_gas_file', 'water_vapor_transmittance_file']

        for key in keys:
            self.l2_prod.attrs[key] = str(self.prod.__dict__[key])
        # rename dim for coarse resolution data
        cams_raster = self.cams.raster.rename({"x": "xc", "y": "yc"})
        transmittance_raster = self.gas_trans.Tg_tot_coarse.rename({"x": "xc", "y": "yc", 'wl': 'wl_all'})

        self.l2_prod['vza'] = native_raster.vza.mean('wl')
        self.l2_prod['sza'] = native_raster.sza
        self.l2_prod['raa'] = native_raster.raa.mean('wl')

        # add cirrus and O2 bands
        cirrus = native_raster.bands.sel(wl=self.wl_cirrus_band)
        cirrus.name = 'cirrus_band'
        cirrus.attrs = {
            'description': 'TOA reflectance in cirrus band from L1C image, might be used to filter out high clouds',
            'wavelength': str(self.wl_cirrus_band)}

        o2band = native_raster.bands.sel(wl=self.wl_o2_band)
        o2band.name = 'o2_band'
        o2band.attrs = {
            'description': 'TOA reflectance in O2 band from L1C image',
            'wavelength': str(self.wl_o2_band)}
        ndwi = native_raster.ndwi
        ndwi_swir = native_raster.ndwi_swir

        # final merge
        self.l2_prod = xr.merge([self.l2_prod,  o2band, cirrus, ndwi, ndwi_swir,
                                 transmittance_raster,cams_raster])

    def to_netcdf(self, ofile):
        '''
        Create output product dimensions, variables, attributes, flags....
        :return:
        '''
        logging.info('export into encoded netcdf')

        complevel = self.complevel

        encoding = {
            'Rrs': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                    "complevel": complevel},
            'Rrs_g': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                      "complevel": complevel},
            'o2_band': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                        "complevel": complevel},
            'cirrus_band': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768,
                            "zlib": True,
                            "complevel": complevel},
            'ndwi': {'dtype': 'int16', 'scale_factor': 0.0001, '_FillValue': -32768, "zlib": True,
                     "complevel": complevel},
            'ndwi_swir': {'dtype': 'int16', 'scale_factor': 0.0001, '_FillValue': -32768, "zlib": True,
                          "complevel": complevel},
            'aot550': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, "zlib": True,
                       "complevel": complevel},
            'vza': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, "zlib": True, "complevel": complevel},
            'raa': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, "zlib": True, "complevel": complevel},
            'sza': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, "zlib": True, "complevel": complevel},
            'aot550': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, "zlib": True,
                       "complevel": complevel},
            'BRDFg': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                      "complevel": complevel},
        }

        # write file
        if os.path.exists(ofile):
            os.remove(ofile)

        odir = os.path.dirname(ofile)
        if not os.path.exists(odir):
            os.mkdir(odir)

        self.l2_prod.to_netcdf(ofile, encoding=encoding)
        self.l2_prod.close()
