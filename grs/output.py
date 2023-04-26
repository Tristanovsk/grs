import os
import xarray as xr
import rioxarray as rio


class l2a_product():
    def __init__(self, prod, l2_prod, cams, gas_trans):
        self.prod = prod
        self.l2_prod = l2_prod
        self.cams = cams
        self.gas_trans = gas_trans

        self.complevel = 5

    def to_netcdf(self,ofile):
        '''
        Create output product dimensions, variables, attributes, flags....
        :return:
        '''

        complevel = self.complevel
        native_raster = self.prod.raster
        # first get attributes from native product
        self.l2_prod.attrs = native_raster.attrs

        keys = ['sunglint_threshold', 'ndwi_threshold', 'green_swir_index_threshold', 'hcld_threshold',
                'dirdata', 'abs_gas_file', 'lut_file', 'water_vapor_transmittance_file']

        for key in keys:
            self.l2_prod.attrs[key] = str(self.prod.__dict__[key])
        # rename dim for coarse resolution data
        cams_raster = self.cams.raster.rename({"x": "xc", "y": "yc"})
        transmittance_raster = self.gas_trans.Tg_tot_coarse.rename({"x": "xc", "y": "yc", 'wl': 'wl_all'}).drop_vars(
            'pressure')
        self.l2_prod['vza'] = native_raster.vza.mean('wl')
        self.l2_prod['sza'] = native_raster.sza
        self.l2_prod['raa'] = native_raster.raa.mean('wl')
        # final merge
        self.l2_prod = xr.merge([self.l2_prod, cams_raster, transmittance_raster])

        encoding = {
            'Rrs': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                    "complevel": complevel},
            'Rrs_g': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
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
        self.l2_prod.to_netcdf(ofile, encoding=encoding)
        self.l2_prod.close()
