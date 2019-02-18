import time,sys
import numpy as np
from esasnappy import ProductIO, GPF, jpy
import geopandas as gpd

#file=sys.argv[1]
file='/home/harmel/Dropbox/rus/S2/gernez/L1C/S2A_MSIL1C_20170210T082051_N0204_R121_T33HYD_20170210T083752.SAFE/'
file='/DATA/Satellite/LANDSAT/L8/L1TP/LC08_L1TP_196029_20180914_20180928_01_T1/LC08_L1TP_196029_20180914_20180928_01_T1_MTL.txt'
shapefile = '/DATA/OBS2CO/vrac/shape/SPO04/SPO04.shp'


def generic_resampler(product, resolution=20, method='Bilinear'):
        '''method: Nearest, Bilinear'''

        GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

        HashMap = jpy.get_type('java.util.HashMap')

        parameters = HashMap()
        parameters.put('targetResolution', resolution)
        parameters.put('upsampling', method)
        parameters.put('downsampling', 'Mean')
        parameters.put('flagDownsampling', 'First')
        parameters.put('resampleOnPyramidLevels', False)

        return GPF.createProduct('Resample', parameters, product)

def get_subset(product, wkt):
    '''subset from wkt POLYGON '''
    SubsetOp = jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
    WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')

    grid = WKTReader().read(wkt)
    op = SubsetOp()
    op.setSourceProduct(product)
    op.setGeoRegion(grid)
    op.setCopyMetadata(True)
    return op.getTargetProduct()




def shp2wkt(shapefile):
    tmp = gpd.GeoDataFrame.from_file(shapefile)
    tmp.to_crs(epsg=4326, inplace=True)
    return  tmp.geometry.values[0].to_wkt()


product = ProductIO.readProduct(file)

p_core = generic_resampler(l2h.product_)

wkt = shp2wkt(shapefile)
p_ = get_subset(p_core,wkt)
w, h = p_.getSceneRasterWidth(),p_.getSceneRasterHeight()
arr__=p_.getBand('green').readPixels(0,0,w,h,np.zeros((w,h)))

#
# def time_load(p_,angles=True):
#
#     tt=time.time
#     tdiff=[]
#     for band in p_.getBandNames(): #('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'):
#         print(band)
#         arr = np.empty(w)
#         t=[]
#         t.append(tt())
#         p_.getBand(band).readPixels(0,0,w,1,arr)
#         #t.append(tt())
#         if angles:
#             angband='view_zenith_'+band
#             p_.getBand(angband).readPixels(0,0,w,1,arr)
#         t.append(tt())
#         tdiff.append(np.diff(t))
#
#     return 'mean time to load a row:  %.4f'%np.mean(tdiff)+'s'
#
# print('1st call with angles from core sampler: ',time_load(p_core))
# print('2nd call with angles from core sampler: ',time_load(p_core))
# print('1st call with angles from s2tbx sampler: ',time_load(p_s2tbx))
# print('2nd call with angles from s2tbx sampler: ',time_load(p_s2tbx))
# print('1st call without angles from s2tbx sampler: ',time_load(p_s2tbx,angles=False))
# print('2nd call without angles from s2tbx sampler: ',time_load(p_s2tbx,angles=False))
#
#
