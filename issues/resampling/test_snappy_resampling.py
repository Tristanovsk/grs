import time,sys
import tracemalloc
import numpy as np
from esasnappy import ProductIO, GPF, jpy


#file=sys.argv[1]
file='/home/harmel/Dropbox/rus/S2/gernez/L1C/S2A_MSIL1C_20170210T082051_N0204_R121_T33HYD_20170210T083752.SAFE/'
file='/DATA/Satellite/LANDSAT/L8/L1TP/LC08_L1TP_196029_20180914_20180928_01_T1/LC08_L1TP_196029_20180914_20180928_01_T1_MTL.txt'


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



product = ProductIO.readProduct(file)
p_core = generic_resampler(product)
p_s2tbx = s2_resampler(product)
w, h = product.getSceneRasterWidth(),product.getSceneRasterHeight()

h=512
def time_load(p_,angles=True):

    tt=time.time
    tdiff=[]
    for band in p_.getBandNames(): #('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'):
        #print(band)
        arr = np.empty((w,h))
        t=[]
        t.append(tt())
        p_.getBand(band).readPixels(0,0,w,h,arr)
        #t.append(tt())
        if angles:
            angband='view_zenith_'+band
            p_.getBand(angband).readPixels(0,0,w,h,arr)
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

# snapshot = tracemalloc.take_snapshot()
# top_stats = snapshot.statistics('lineno')
#
# print("[ Top 10 ]")
# for stat in top_stats[:10]:
#     print('s2tbx',stat)
# print('1st call with angles from core sampler: ',time_load(p_core))
# print('2nd call with angles from core sampler: ',time_load(p_core))
# print('1st call with angles from s2tbx sampler: ',time_load(p_s2tbx))
# print('2nd call with angles from s2tbx sampler: ',time_load(p_s2tbx))
