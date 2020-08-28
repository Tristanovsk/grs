import time,sys
import numpy as np
from esasnappy import ProductIO, GPF, jpy

from grs import s2_angle

file='/DATA/Satellite/SENTINEL2/venice/L1C/S2A_OPER_PRD_MSIL1C_PDMC_20170121T221023_R022_V20150902T101026_20150902T101026.SAFE'

sampler=sys.argv[1]
file=sys.argv[2]


def generic_resampler(product, resolution=20, method='Bilinear'):
        '''method: Nearest, Bilinear'''

        GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

        HashMap = jpy.get_type('java.util.HashMap')

        parameters = HashMap()
        parameters.put('targetResolution', resolution)
        parameters.put('upsampling', method)
        parameters.put('downsampling', 'Mean')
        parameters.put('flagDownsampling', 'FlagMedianAnd')
        parameters.put('resampleOnPyramidLevels', False)

        return GPF.createProduct('Resample', parameters, product)

def s2_resampler(product, resolution=20, upmethod='Bilinear', downmethod='Mean',
                     flag='FlagMedianAnd', opt=False):

        '''Resampler operator dedicated to Sentinel2-msi characteristics (e.g., viewing angles)
        :param product: S2-msi product as provided by esasnappy ProductIO.readProduct()
        :param resolution: target resolution in meters (10, 20, 60)
        :param upmethod: interpolation method ('Nearest', 'Bilinear', 'Bicubic')
        :param downmethod: aggregation method ('First', 'Min', 'Max', 'Mean', 'Median')
        :param flag: method for flags aggregation ('First', 'FlagAnd', 'FlagOr', 'FlagMedianAnd', 'FlagMedianOr')
        :param opt: resample on pyramid levels (True/False)
        :return: interpolated target product'''

        res = str(resolution)

        resampler = jpy.get_type('org.esa.s2tbx.s2msi.resampler.S2ResamplingOp')

        op = resampler()
        op.setParameter('targetResolution', res)
        op.setParameter('upsampling', upmethod)
        op.setParameter('downsampling', downmethod)
        op.setParameter('flagDownsampling', flag)
        op.setParameter('resampleOnPyramidLevels', opt)
        op.setSourceProduct(product)

        return op.getTargetProduct()

#------------
# calculation
import os
xml=os.path.join(file,'GRANULE/S2A_OPER_MSI_L1C_TL_EPA__20161008T200750_A001020_T32TQR_N02.04/S2A_OPER_MTD_L1C_TL_EPA__20161008T200750_A001020_T32TQR.xml')



s2_angle.s2angle_v2().angle_computation(xml)

product = ProductIO.readProduct(file)
product.getMetadataRoot().getElement('Granules').getElement('Level-1C_Tile_12QWM').getElement('Geometric_Info')


p_core = generic_resampler(product)
p_s2tbx = s2_resampler(product)
w, h = product.getSceneRasterWidth(),product.getSceneRasterHeight()

h=512
def time_load(p_,angles=True):

    tt=time.time
    tdiff=[]
    for band in p_.getBandNames():
        arr = np.empty((w,h))
        t=[]
        t.append(tt())
        arr = p_.getBand(band).readPixels(0,0,w,h,arr)
        t.append(tt())
        tdiff.append(np.diff(t))

    return 'mean time to load a row:  %.4f'%np.mean(tdiff)+'s'

arg=sys.argv[1]
if arg == 'core':
    print('1st call with angles from core sampler: ',time_load(p_core,angles=False))
    print('2nd call with angles from core sampler: ',time_load(p_core,angles=False))
else:
    print('1st call without angles from s2tbx sampler: ',time_load(p_s2tbx,angles=False))
    print('2nd call without angles from s2tbx sampler: ',time_load(p_s2tbx,angles=False))
