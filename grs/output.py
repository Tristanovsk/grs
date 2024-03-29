'''
Module to set output product and format.
'''
import os
import numpy as np
import xarray as xr
import logging


class L2aProduct():
    '''
    Create and handle L2A object
    '''
    def __init__(self, prod, l2_prod, cams, gas_trans,dem=None):
        self.prod = prod
        self.l2_prod = l2_prod
        self.cams = cams
        self.gas_trans = gas_trans
        self.dem = dem
        self.ancillary = None
        self.wl_cirrus_band = 1375
        self.wl_wv_band = 945
        self.complevel = 5
        self.construct_l2a()

    def construct_l2a(self):
        '''
        Construct rioxarray from raster data and metadata as attributes.

        :return:
        '''
        logging.info('construct l2a')
        native_raster = self.prod.raster
        # first get attributes from native product
        self.l2_prod.attrs = native_raster.attrs

        keys = ['sunglint_threshold', 'ndwi_threshold', 'vis_swir_index_threshold', 'hcld_threshold',
                'dirdata', 'abs_gas_file', 'water_vapor_transmittance_file']

        for key in keys:
            self.l2_prod.attrs[key] = str(self.prod.__dict__[key])
        # rename dim for coarse resolution data
        cams_raster = self.cams.raster.rename({"x": "xc", "y": "yc"})
        transmittance_raster = self.gas_trans.Tg_tot_coarse.rename({"x": "xc", "y": "yc", 'wl': 'wl_all'})

        # save geometrie angles, take mean values if angle depends on spectral bands
        def angle_extraction(angle):
            if 'wl' in angle.coords:
                return angle.mean('wl')
            else:
                return angle
        self.l2_prod['vza'] = angle_extraction(native_raster.vza)
        self.l2_prod['sza'] = native_raster.sza
        self.l2_prod['raa'] = angle_extraction(native_raster.raa)

        # add cirrus and water vapor bands
        if self.prod.bcirrus:
            cirrus = self.prod.cirrus
            cirrus.name = 'cirrus_band'
            cirrus.attrs = {
                'description': 'TOA reflectance in cirrus band from L1C image, might be used to filter out high clouds',
                'wavelength': str(self.wl_cirrus_band)}

        if self.prod.bwv:
            wvband = self.prod.wv
            wvband.name = 'wv_band'
            wvband.attrs = {
                'description': 'TOA reflectance in water vapor band from L1C image',
                'wavelength': str(self.wl_wv_band)}

        ndwi = native_raster.ndwi
        ndwi_swir = native_raster.ndwi_swir

        # ------------------
        # adding digital elevation model to output
        dem =self.dem
        # add empty DEM in raster
        if dem is None:
            dem = xr.full_like(self.l2_prod['sza'], np.nan, dtype=np.float32)
            dem.name = 'dem'

        # final merge
        self.l2_prod = xr.merge([self.l2_prod, native_raster.flags_l1c,dem,native_raster.surfwater])
        self.l2_prod.rio.set_spatial_dims(x_dim='x', y_dim='y', inplace=True)
        self.l2_prod.rio.write_coordinate_system(inplace=True)
        self.l2_prod.rio.write_crs(inplace=True)

        self.ancillary = xr.merge([transmittance_raster, cams_raster])
        self.ancillary.rio.set_spatial_dims(x_dim='xc', y_dim='yc', inplace=True)
        self.ancillary.rio.write_coordinate_system(inplace=True)
        self.ancillary.rio.write_crs(inplace=True)

    def to_netcdf(self, output_path, snap_compliant=False):
        '''
        Create output product dimensions, variables, attributes, flags....

        :return:
        '''
        logging.info('export into encoded netcdf')

        complevel = self.complevel

        encoding = {
            'aot550': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999, "zlib": True,
                       "complevel": complevel, 'grid_mapping': 'spatial_ref'},
            'BRDFg': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                      "complevel": complevel, 'grid_mapping': 'spatial_ref'},

            'dem': {'dtype': 'int16', 'scale_factor': 0.00001,'add_offset': 1000, '_FillValue': -32768, "zlib": True,
                     "complevel": complevel, 'grid_mapping': 'spatial_ref'},

            'surfwater': {'dtype': 'int8',  "complevel": complevel, "zlib": True, 'grid_mapping': 'spatial_ref'},
            'vza': {'dtype': 'int16', 'scale_factor': 0.01, '_FillValue': -32768, "zlib": True, "complevel": complevel,
                    'grid_mapping': 'spatial_ref'},
            'raa': {'dtype': 'int16', 'scale_factor': 0.01, 'add_offset': 180, '_FillValue': -32768, "zlib": True,
                    "complevel": complevel, 'grid_mapping': 'spatial_ref'},
            'sza': {'dtype': 'int16', 'scale_factor': 0.01, 'add_offset': 30, '_FillValue': -32768,
                    "complevel": complevel, 'grid_mapping': 'spatial_ref'},
            'flags_l1c': {"complevel": complevel, 'grid_mapping': 'spatial_ref'}
        }

        if snap_compliant:
            # explode Rrs array to get one variable per band to be SNAP compliant
            for var in ['Rrs']:
                img_snap = self.l2_prod[var].to_dataset("wl")
                suff = ''
                if var == 'Rrs_g':
                    suff = 'with sunglint '
                for ii,band in enumerate(img_snap.keys()):
                    band_name = var + '_{:d}'.format(band)
                    band_ref = 'B{:d}'.format(band)
                    wavelength = self.prod.raster.wl_true.values[ii]
                    #bandwidth = self.prod.l1c.band_info[band_ref]['bandwidth']
                    img_snap[band].attrs = {
                        'long_name': 'Remote sensing reflectance' + suff + ' at {:d} nm'.format(band),
                        'Unit': 'sr-1',
                        'units': 'sr-1',
                        'radiation_wavelength': wavelength,
                        'radiation_wavelength_unit': 'nm',
                        #'bandwidth': bandwidth,
                        'wavelength': wavelength,
                        'valid_pixel_expression': ''}
                    encoding[band_name] = {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3,
                                           '_FillValue': -32768, "zlib": True,
                                           "complevel": complevel, 'grid_mapping': 'spatial_ref'}

                    img_snap = img_snap.rename({band: band_name})
                self.l2_prod = xr.merge([self.l2_prod.drop_vars(var), img_snap])
            self.l2_prod.attrs['auto_grouping'] = 'Rrs:Rrs_g'
            self.l2_prod.attrs['metadata_profile'] = 'beam'
        else:
            for band_name in ['Rrs']:
                encoding[band_name] = {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3,
                                       '_FillValue': -32768, "zlib": True,
                                       "complevel": complevel, 'grid_mapping': 'spatial_ref'}


            self.l2_prod['central_wavelength'] = ('wl', self.prod.raster.wl_true.values)
            self.l2_prod = self.l2_prod.set_coords('central_wavelength')
            self.l2_prod.attrs['metadata_profile'] = 'datacube'

        # get file naming and create folder
        basename = os.path.basename(output_path)
        ofile = os.path.join(output_path, basename)

        # create directory if not existing
        os.makedirs(output_path,exist_ok=True)

        # clean up to avoid permission denied
        if os.path.exists(ofile + '.nc'):
            os.remove(ofile + '.nc')
        if os.path.exists(ofile + '_anc.nc'):
            os.remove(ofile + '_anc.nc')

        # export full raster data
        # TODO check why or generalize the following approach:
        # fix for conflicts with attrs and encoding, needs to remove 'grid_mapping' from input attrs
        self.l2_prod.sza.attrs=''
        self.l2_prod.surfwater.attrs=''

        self.l2_prod.to_netcdf(ofile + '.nc', encoding=encoding)

        # self.l2_prod.close()

        # export ancillary data (coarse resolution)
        encoding = {}
        for variable in list(self.ancillary.keys()):
            encoding[variable] = {"zlib": True, "complevel": complevel, 'grid_mapping': 'spatial_ref'}

        self.ancillary.to_netcdf(ofile + '_anc.nc', 'w', encoding=encoding)  # ,group='ancillary')

        # self.ancillary.close()

        return
