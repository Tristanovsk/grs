'''
Module dedicated to pixel classification.
'''

import numpy as np
import xarray as xr

import logging

from s2cloudless import S2PixelCloudDetector

FLAG_NAME = 'flags_l1c'

class Settings:
    '''
    Settings used for the different masking procedures.
    '''
    def __init__(self):
        self.bvis = 490
        self.bnir = 842
        self.bswir = 1610
        self.bswir2 = 2190

        # set extra bands (e.g., cirrus, water vapor)
        self.bcirrus = 1375
        self.bwv = 945

        # mask thresholding parameters
        self.swir_threshold = 0.2
        self.ndwi_threshold = 0.0
        self.vis_swir_index_threshold = 0.
        self.thin_cirrus_threshold = 3e-3
        self.opac_cirrus_threshold = 9e-3

        # settings for s2cloudless masks
        self.low_confid_cloud_proba_thresh = 0.6
        self.low_confid_cloud_dilation = 8
        self.high_confid_cloud_proba_thresh = 0.8
        self.high_confid_cloud_dilation = 5

        self.number_of_flags = 32


class Masking(Settings):
    '''
    Class dedicated to pixel classification and masking.
    '''
    def __init__(self, prod):
        '''

        :param prod: image raster containing the spectral bands
        '''
        Settings.__init__(self)
        self.prod = prod
        self.flag_descriptions = np.empty(self.number_of_flags,
                                          dtype='object')#np.dtype('U', 1000))
        self.flag_names = np.empty(self.number_of_flags,
                                   dtype='object')#np.dtype('U', 1000))

        # generate flag raster from nodata mask
        self.create_flags()
        self.flags.attrs['long_name'] = 'flags computed from l1c image'
        self.flag_stats = {}

    def create_flags(self,
                     bitmask=0,
                     name='nodata',
                     description='nodata in input image'):
        '''
        set flag nodata: condition either nan or crazy numerical values

        :return:
        '''

        self.flags = (np.isnan(self.prod.bands.isel(wl=1)) | (self.prod.bands.isel(wl=1) > 1e3)).astype(np.int32)
        # set naming and description as attributes
        self.flag_descriptions[bitmask] = description
        self.flag_names[bitmask] = name

    def cloud_mask(self,
                   bitmasks=[1, 2],
                   names=['cloud_p06', 'cloud_p08'],
                   descriptions=['low confidence cloud mask from s2cloudless with settings proba.',
                                 'high confidence cloud mask from s2cloudless with settings proba.']
                   ):
        '''
        Apply s2cloudless masking with two levels of confidence

        :param bitmasks: bit number on which the flag is coded
        :param names: name of the flag
        :param descriptions: description of the flag
        :return:
        '''

        logging.info('cloud masking with s2cloudless')

        # apply cloud mask
        bands = np.array([self.prod.bands.values.transpose((1, 2, 0))])
        # ---------------------
        # first low confidence
        # ---------------------
        bitmask = bitmasks[0]

        cloud_detector = S2PixelCloudDetector(threshold=self.low_confid_cloud_proba_thresh, average_over=1,
                                              dilation_size=self.low_confid_cloud_dilation, all_bands=True)
        probability_maps = cloud_detector.get_cloud_probability_maps(bands)
        cloud_mask = cloud_detector.get_mask_from_prob(probability_maps)[0]
        self.flags = self.flags + ((cloud_mask == 1) << bitmask)

        # add name and description
        self.flag_descriptions[bitmask] = descriptions[0] + ' threshold={:.2f}, dilation_size={:d}'.format(
            self.low_confid_cloud_proba_thresh, self.low_confid_cloud_dilation)
        self.flag_names[bitmask] = names[0]

        # ---------------------
        # second high confidence
        # ---------------------
        bitmask = bitmasks[1]

        cloud_detector = S2PixelCloudDetector(threshold=self.high_confid_cloud_proba_thresh, average_over=1,
                                              dilation_size=self.high_confid_cloud_dilation, all_bands=True)
        cloud_mask = cloud_detector.get_mask_from_prob(probability_maps)[0]
        self.flags = self.flags + ((cloud_mask == 1) << bitmask)

        # add name and description
        self.flag_descriptions[bitmask] = descriptions[1] + ' threshold={:.2f}, dilation_size={:d}'.format(
            self.high_confid_cloud_proba_thresh, self.high_confid_cloud_dilation)
        self.flag_names[bitmask] = names[1]

        del cloud_detector, cloud_mask, probability_maps

    def land_mask(self,
                  bitmasks=[3, 4],
                  names=['land_swir_visible_index', 'land_red_visible_index'],
                  descriptions=['land-water mask from normalized index from swir and visible band',
                                'land-water mask from normalized index from nir and visible band, warning could fail for turbid waters'],
                  ):
        '''
        apply water/land mask
        compute mask from NDWI from visible and NIR and visible and SWIR

        :param bitmasks: bit numbers on which the flags are coded
        :param names: name of the flags
        :param descriptions: description of the flags
        :return:
        '''

        logging.info('land masking')
        visible = self.prod.bands.sel(wl=self.bvis, method='nearest')
        nir = self.prod.bands.sel(wl=self.bnir, method='nearest')
        swir = self.prod.bands.sel(wl=self.bswir, method='nearest')

        ndwi = (visible - nir) / (visible + nir)
        ndwi_swir = (visible - swir) / (visible + swir)

        # set flags raster
        self.flags = self.flags + (
                ((ndwi_swir.values < self.vis_swir_index_threshold) << bitmasks[0]) +
                ((ndwi.values < self.ndwi_threshold) << bitmasks[1]))

        # set naming and description as attributes
        self.flag_descriptions[bitmasks[0]] = (descriptions[0] +
                                                ', bands centered on {:.1f} and {:.1f} nm'.format(self.bvis,self.bswir)+
                                               ', threshold={:.2f}'.format(self.vis_swir_index_threshold))
        self.flag_names[bitmasks[0]] = names[0]

        self.flag_descriptions[bitmasks[1]] = (descriptions[1] +
                                               ', bands centered on {:.1f} and {:.1f} nm'.format(self.bvis, self.bnir) +
                                               ', threshold={:.2f}'.format(self.ndwi_threshold))
        self.flag_names[bitmasks[1]] = names[1]

    def high_swir_mask(self,
                       bitmask=7,
                       name='high_swir',
                       description='high swir (bright cloud, too bright reflection...)'
                       ):
        '''
        Flag to mask (too) high values of the swir bands, generally to remove potential cloud cotamination.

        :param bitmask: bit number on which the flags are coded
        :param name: name of the flag
        :param description: description used as attribute
        :return:
        '''

        logging.info('high swir masking')
        b2200 = self.prod.bands.sel(wl=self.bswir2)
        self.flags = self.flags + ((b2200.values > self.swir_threshold) << bitmask)
        # set naming and description as attributes
        self.flag_descriptions[bitmask] = description +', threshold {:.4f}'.format(self.swir_threshold)

        self.flag_names[bitmask] = name

    def cirrus_mask(self,
                    bitmasks=[5,6],
                    names=['thin_cirrus','opac_cirrus'],
                    descriptions=['thin cirrus mask from cirrus band',
                                  'opac cirrus mask from cirrus band']
                    ):
        '''
        Compute cirrus mask for thin and opac high clouds from spectral band around 1375 nm, if exists

        :param bitmasks: bit numbers on which the flags are coded
        :param names: name of the flags
        :param descriptions: description of the flags
        :return:
        '''

        logging.info('cirrus masking')
        cirrus = self.prod.bands.sel(wl=self.bcirrus)
        thresholds = [self.thin_cirrus_threshold,self.opac_cirrus_threshold]
        for ii in range(2):
            self.flags = self.flags + ((cirrus.values > thresholds[ii]) << bitmasks[ii])

            # set naming and description as attributes
            self.flag_descriptions[bitmasks[ii]] = descriptions[ii] +', threshold {:.4f}'.format(thresholds[ii])
            self.flag_names[bitmasks[ii]] = names[ii]

    def get_stats(self):
        '''
        Compute image statistics for each flag and save them into dictionary

        :return:
        '''

        flag_value = 1
        for ii, flag_name in enumerate(self.prod[FLAG_NAME].flag_names):
            if flag_name == 'None':
                continue
            flag = ((self.prod[FLAG_NAME] & flag_value) != 0)
            flag_stat = float(flag.sum() / flag.count())
            self.flag_stats['flag_' + flag_name] = flag_stat
            flag_value = flag_value << 1


    def process(self,
                output="prod"):
        '''
        Generate the flags raster and attributes

        :param output:
            - if None returns nothing but the raster is updated within the masking object
            - if "prod" returns the full raster updated with the flags variable and attributes
            - if "flags" returns xarray DataArray of the flags
              plus the dictionary of flags statistics
        :return: see :param output
        '''
        # apply the masking processors
        self.cloud_mask()
        self.land_mask()
        self.cirrus_mask()
        self.high_swir_mask()

        if isinstance(self.flags, xr.DataArray):
            self.prod[FLAG_NAME] = (('y', 'x'), self.flags.values)
        else:
            self.prod[FLAG_NAME] = (('y', 'x'), self.flags)

        # set the attributes with the names and description of the flags
        self.prod[FLAG_NAME].attrs['flag_descriptions'] = self.flag_descriptions.astype(str)
        self.prod[FLAG_NAME].attrs['flag_names'] = self.flag_names.astype(str)

        # compute flag statistics over the image
        self.get_stats()
        self.prod.attrs.update(self.flag_stats)

        if output == "prod":
            return self.prod
        elif output == "flags":
            return self.prod[FLAG_NAME], self.flag_stats
