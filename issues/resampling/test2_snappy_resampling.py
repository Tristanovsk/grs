import time,sys
import numpy as np
from esasnappy import ProductIO, GPF, jpy


sampler=sys.argv[1]
file=sys.argv[2]

def generic_resampler(product, resolution=20, method='Nearest'):
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

def s2_resampler(product, resolution=20, upmethod='Nearest', downmethod='Mean',
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


file='/DATA/Satellite/SENTINEL2/SPO04/L1C/S2B_MSIL1C_20180927T103019_N0206_R108_T31TGK_20180927T143835.SAFE'
product = ProductIO.readProduct(file)
bands = product.getBandNames()
p_core = generic_resampler(product)
p_s2tbx = s2_resampler(product)
w, h = product.getSceneRasterWidth(),product.getSceneRasterHeight()

h=512
def time_load(p_,angles=True):

    tt=time.time
    tdiff=[]
    for band in list(bands)[::-1]:
        arr = np.empty((w,h))
        t=[]
        t.append(tt())
        p_.getBand(band).readPixels(0,0,w,h,arr)
        p_.getBand(band).dispose()
        t.append(tt())
        tdiff.append(np.diff(t))

    return 'mean time to load a row:  %.4f'%np.mean(tdiff)+'s'


if sampler == 'core':
    print('1st call with angles from core sampler: ',time_load(p_core,angles=False))
    print('2nd call with angles from core sampler: ',time_load(p_core,angles=False))
else:
    print('1st call without angles from s2tbx sampler: ',time_load(p_s2tbx,angles=False))
    print('2nd call without angles from s2tbx sampler: ',time_load(p_s2tbx,angles=False))
