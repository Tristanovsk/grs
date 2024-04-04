import os
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib

from grs import AuxData

lat, lon = 47, -1.5

cams_folder = os.path.abspath('/nfs/DP/ECMWF/CAMS')
type='cams_forecast'
#type='cams_reanalysis'

def datespan(startDate, endDate, delta=relativedelta(months=+1)):
    currentDate = startDate
    while currentDate < endDate:
        yield currentDate
        currentDate += delta

for date in datespan(dt.datetime(2019,7,17,12),dt.datetime.now()+relativedelta(months=-1)):#.datetime(2019,8,17,12)):
    print(date)
    cams_file=os.path.join(cams_folder, date.strftime('%Y-%m') + '_month_'+type+'.nc')

    cams = AuxData.cams()
    cams.load_cams_data(cams_file, date, data_type=type)

def wktbox(center_lon, center_lat, width=1, height=1):
    '''

    :param center_lon: decimal longitude
    :param center_laat: decimal latitude
    :param width: width of the box in km
    :param height: haight of the box in km
    :return: wkt of the box centered on provided coordinates
    '''
    from math import sqrt, atan, pi
    import pyproj
    geod = pyproj.Geod(ellps='WGS84')
    width, height = width * 1000, height * 1000
    rect_diag = sqrt(width ** 2 + height ** 2)

    azimuth1 = atan(width / height)
    azimuth2 = atan(-width / height)
    azimuth3 = atan(width / height) + pi  # first point + 180 degrees
    azimuth4 = atan(-width / height) + pi  # second point + 180 degrees

    pt1_lon, pt1_lat, _ = geod.fwd(center_lon, center_lat, azimuth1 * 180 / pi, rect_diag)
    pt2_lon, pt2_lat, _ = geod.fwd(center_lon, center_lat, azimuth2 * 180 / pi, rect_diag)
    pt3_lon, pt3_lat, _ = geod.fwd(center_lon, center_lat, azimuth3 * 180 / pi, rect_diag)
    pt4_lon, pt4_lat, _ = geod.fwd(center_lon, center_lat, azimuth4 * 180 / pi, rect_diag)

    wkt_point = 'POINT (%.6f %.6f)' % (center_lon, center_lat)
    wkt_poly = 'POLYGON (( %.6f %.6f, %.6f %.6f, %.6f %.6f, %.6f %.6f, %.6f %.6f ))' % (
        pt1_lon, pt1_lat, pt2_lon, pt2_lat, pt3_lon, pt3_lat, pt4_lon, pt4_lat, pt1_lon, pt1_lat)
    return wkt_poly

wkt = wktbox(lon, lat,50,50)


cams.get_cams_ancillary(cams_file, date, wkt)
for i in range(4):
    date=date + datetime.timedelta(days=1)
    cams.get_cams_aerosol(cams_file, date, wkt)
    print(date.ctime(),cams.aot_rast[1,...].min(),cams.aot_rast[1,...].mean(),cams.aot_rast[1,...].max())


