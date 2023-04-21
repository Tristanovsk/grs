'''
Main program
'''

from pathlib import Path
import os, shutil
import zipfile
import tarfile
import glob

import matplotlib

matplotlib.use('tkagg')
import matplotlib.pyplot as plt

plt.figure()
import numpy as np
import xarray as xr
import logging

from s2driver import driver_S2_SAFE as S2
from grs import product, acutils, utils

# from grs import product, l2a, utils, acutils, auxdata

# from grs.fortran.grs import main_algo as grs_solver

file = '/sat_data/satellite/sentinel2/L1C/31TFJ/S2A_MSIL1C_20201004T104031_N0209_R008_T31TFJ_20201004T125253.SAFE'
file = '/media/harmel/vol1/Dropbox/satellite/S2/L1C/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE'
file_nc = file.replace('.SAFE', '.nc')
bandIds = range(13)
resolution = 60

if not os.path.exists(file_nc):
    l1c = S2.s2image(file, band_idx=bandIds, resolution=resolution)
    print(l1c.crs)
    l1c.load_product()
    prod = product(l1c.prod)
    encoding = {'bands': {'dtype': 'int16', 'scale_factor': 0.00001,'add_offset':.3, '_FillValue': -32768},
                'vza': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999},
                'raa': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999},
                'sza': {'dtype': 'int16', 'scale_factor': 0.001, '_FillValue': -9999}}
    l1c.prod.to_netcdf(file_nc, encoding=encoding)
    l1c.prod.close()

else:
    prod = product(xr.open_dataset(file_nc))

##################################
# Fetch optional mask products
# resample for common resolution
# subset to ROI
##################################
# TODO activate MAJA reader and flag retrieval
# prod.get_flags()
# logging.info('fetching flags...')
# maja, waterdetect = None, None
# if maja_xml:
#     try:
#         maja = ProductIO.readProduct(maja_xml)
#         maja = _utils.resampler(maja, resolution=resolution)
#         maja = _utils.get_subset(maja, wkt)
#
#     except:
#         logging.info('!!! issues with ' + maja_xml + '; please check if file exists')
#         raise
##################################
# GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
##################################
# prod.get_cams()
##################################
## ADD ELEVATION AND PRESSURE BAND
##################################
# prod.get_elevation()
# logging.info('adding elevation band...')
# if dem:
#     logging.info('add elevation band')
#     # ~ high_latitude = (latmax >= 60) | (latmin <= -60)
#     prod.get_elevation(dem_file=dem_file, dem_glo30_dir=dem_glo30_dir)
# else:
#     prod.elevation = np.zeros([prod.height, prod.width])

#####################################
# LOAD LUT FOR ATMOSPHERIC CORRECTION
#####################################
logging.info('loading lut...' + prod.lutfine)
lutf = acutils.lut(prod.band_names)
lutc = acutils.lut(prod.band_names)
lutf.load_lut(prod.lutfine, prod.sensordata.indband)
lutc.load_lut(prod.lutcoarse, prod.sensordata.indband)

# reproject lut array on the angles of the image
# angles are rounded to reduce the dims of interpolated LUT
_utils = utils()
sza_ = _utils.remove_na(np.unique(prod.raster.sza.round(1)))
vza_ = _utils.remove_na(np.unique(prod.raster.vza.round(1)))
azi_ = _utils.remove_na(np.unique(prod.raster.raa.round(0)))
lutf.interp_n_slice(sza_, vza_, azi_)
lutc.interp_n_slice(sza_, vza_, azi_)
aotlut = np.array(lutf.aot, dtype=prod.type)

####################################
#    absorbing gases correction
####################################
gas_trans=acutils.gaseous_transmittance(prod,prod.gas_lut,prod.Twv_lut)
gas_trans.get_gaseous_transmittance()

logging.info('loading SMAC algorithm...')
smac = acutils.smac(prod.sensordata.smac_bands, prod.sensordata.smac_dir)
smac.set_gas_param()
smac.set_values(o3du=prod.cams.o3du, h2o=prod.cams.h2o)
smac.set_standard_values(prod.cams.pressure_msl)
prod.cams.no2 = smac.uno2


######################################
# Water mask
######################################
# Compute NDWI
green = prod.raster.bands.sel(wl=prod.b565)
nir = prod.raster.bands.sel(wl=prod.b865)
swir = prod.raster.bands.sel(wl=prod.b1600)
b2200 = prod.raster.bands.sel(wl=prod.b2200)

ndwi = (green - nir) / (green + nir)
ndwi_swir = (green - swir) / (green + swir)
self=prod
masked_raster = prod.raster.bands.where(ndwi > self.ndwi_threshold). \
            where(b2200 < self.sunglint_threshold). \
            where(ndwi_swir > self.green_swir_index_threshold)





from matplotlib.colors import ListedColormap
bcmap = ListedColormap(['khaki', 'lightblue'])


def water_mask(ndwi, threshold=0):
    water = xr.where(ndwi > threshold, 1, 0)
    return water.where(~np.isnan(ndwi))


def plot_water_mask(ndwi, ax, threshold=0):
    water = water_mask(ndwi, threshold)

    water.plot.imshow(ax=ax, cmap=bcmap,
                      cbar_kwargs={'ticks': [0, 1], 'shrink': shrink})
    ax.set_title(str(threshold) + ' < NDWI')


ndwi_ = ndwi
fig, axs = plt.subplots(2, 2, figsize=(17, 15), sharex=True, sharey=True)
fig.subplots_adjust(bottom=0.1, top=0.95, left=0.1, right=0.99,
                    hspace=0.05, wspace=0.05)
shrink = 0.8
axs = axs.ravel()

fig = ndwi_.plot.imshow(ax=axs[0], cmap=plt.cm.BrBG, robust=True,
                        cbar_kwargs={'shrink': shrink})
# axes.coastlines(resolution='10m',linewidth=1)
axs[0].set_title('Sentinel 2, NDWI')

for i, threshold in enumerate([-0.2, 0., 0.2]):
    plot_water_mask(ndwi_, axs[i + 1], threshold=threshold)

plt.show()

lonmin, lonmax, latmin, latmax = 4.98, 5.23, 43.395, 43.56
wkt = "POLYGON((" + str(lonmax) + " " + str(latmax) + "," + str(lonmax) + " " \
      + str(latmin) + "," + str(lonmin) + " " + str(latmin) + "," + str(lonmin) + " " \
      + str(latmax) + "," + str(lonmax) + " " + str(latmax) + "))"

sensor = None

ancillary = 'cds_forecast'
aerosol = ancillary
altitude = 0,
dem = True
aeronet_file = None
aot550 = 0.1
angstrom = 1
resolution = 20
unzip = False
untar = False
allpixels = False
waterdetect_file = None
waterdetect_only = False
memory_safe = False
angleonly = False
grs_a = False
output = 'Rrs'

##################################
# Generate prod object
##################################

_utils = utils()
sensor = _utils.get_sensor(file)
sensordata = auxdata.sensordata(sensor)
resolution = sensordata.resolution
indband = sensordata.indband

prod = info(l1c, sensordata, aerosol, ancillary, output)
prod.headerfile = file

##################################
# Raw water mask
##################################
# Compute NDWI
green = l1c.bands.sel(wl=560)
nir = l1c.bands.sel(wl=865)
ndwi = (green - nir) / (green + nir)
threshold = 0
l1c.bands = l1c.bands.where(ndwi > threshold)

coarsening = 1
plt.figure()
l1c.bands[:, ::coarsening, ::coarsening].plot.imshow(col='wl', col_wrap=3, vmax=0.5, cmap=plt.cm.binary_r,
                                                     subplot_kws=dict(xlabel='', xticks=[], yticks=[]))


def water_mask(ndwi, threshold=0):
    water = xr.where(ndwi > threshold, 1, 0)
    return water.where(~np.isnan(ndwi))


##################################
# GET METADATA
##################################
logging.info('getting metadata...')
# TODO clean up this part and other metadata to be loaded
if 'S2' in sensor:
    meta = prod.product.getMetadataRoot().getElement('Level-1C_User_Product').getElement(
        'General_Info').getElement(
        'Product_Image_Characteristics').getElement('Reflectance_Conversion')
    prod.U = float(str(meta.getAttribute('U').getData()))
    prod.solar_irr = np.zeros(len(indband), dtype=np.float32)
    for i, iband in zip(range(len(indband)), indband):
        prod.solar_irr[i] = float(str(meta.getElement('Solar_Irradiance_List').getAttributeAt(iband).getData()))

else:
    meta = prod.product.getMetadataRoot().getElement("LANDSAT_METADATA_FILE").getElement("IMAGE_ATTRIBUTES")
    prod.U = float(str(meta.getAttribute('EARTH_SUN_DISTANCE').getData())) ** 2
    prod.solar_irr = np.array(prod.sensordata.solar_irr)[indband]

# convert into mW cm-2 um-1
prod.solar_irr = prod.solar_irr / 10

##################################
# SUBSET TO AREA OF INTEREST
##################################


##################################
# Fetch optional mask products
# resample for common resolution
# subset to ROI
##################################
logging.info('fetching flags...')
maja, waterdetect = None, None
# if maja_xml:
#     try:
#         maja = ProductIO.readProduct(maja_xml)
#         maja = _utils.resampler(maja, resolution=resolution)
#         if wkt is not None:
#             maja = _utils.get_subset(maja, wkt)
#
#     except:
#         logging.info('!!! issues with ' + maja_xml + '; please check if file exists')
#         raise
#
# if waterdetect_file:
#     try:
#         waterdetect = ProductIO.readProduct(waterdetect_file)
#         waterdetect = _utils.get_subset(waterdetect, wkt)
#         waterdetect = _utils.resampler(waterdetect, resolution=resolution)
#
#     except:
#         logging.info('!!! issues with ' + waterdetect_file + '; please check if file exists')
#         raise

##################################
## ADD ELEVATION BAND
##################################
logging.info('adding elevation band...')
if dem:
    logging.info('add elevation band')
    high_latitude = (latmax >= 60) | (latmin <= -60)
    prod.get_elevation(high_latitude)

else:
    prod.elevation = np.zeros([prod.height, prod.width])

##################################
# GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
##################################
logging.info('getting CAMS data...')
prod.aux = auxdata.cams()

if ancillary != 'default':
    if prod.aerosol == 'cds_forecast':
        target = os.path.join(prod.cams_folder, prod.date.strftime('%Y'), prod.date.strftime('%Y-%m') +
                              '_month_cams-global-atmospheric-composition-forecasts.nc')
        prod.aux.get_cams_ancillary(target, prod.date, prod.wkt, param=['msl', 'gtco3', 'tcwv', 'tcno2', 't2m'])
    else:
        target = Path(os.path.join(prod.cams_folder, prod.date.strftime('%Y'),
                                   prod.date.strftime('%Y-%m') + '_month_' + prod.ancillary + '.nc'))
        # do not load here since already implemented elsewhere in CNES HPC
        # prod.aux.load_cams_data(target, prod.date, data_type=prod.ancillary)
        prod.aux.get_cams_ancillary(target, prod.date, prod.wkt)

## uncomment this part to use ecmwf files provided in the .SAFE format
# if 'S2' in sensor:
#     prod.aux.get_tile_dir(file)
#     prod.aux.get_aux_dir()
#     prod.aux.get_ecmwf_data()

# get pressure at the scene altitude
prod.pressure_msl = prod.aux.msl  # acutils.misc.get_pressure(altitude, prod.aux.msl)
if dem:
    altitude = prod.elevation
    altitude[altitude < -20] = 0
prod.pressure = acutils.misc.get_pressure(altitude, prod.pressure_msl)  # prod.aux.pressure = prod.pressure

######################################
#      Create output l2 product
#          'l2_product'
######################################
logging.info('creating L2 output product')
prod.create_product(maja=maja, waterdetect=waterdetect)
try:
    prod.load_data()
except:
    if unzip:
        # remove unzipped files (Sentinel files)
        shutil.rmtree(file, ignore_errors=True)
    if untar:
        # remove untared files (Landsat files)
        shutil.rmtree(tmp_dir, ignore_errors=True)
    raise NameError('No data available for requested area')
prod.load_flags()

# ---------
# delete Cache
# TODO need to find a way to delete the processed file only (keep the others)
# prod.deleteCache()

#####################################
# LOAD LUT FOR ATMOSPHERIC CORRECTION
#####################################
logging.info('loading lut...' + prod.lutfine)
lutf = acutils.lut(prod.band_names)
lutc = acutils.lut(prod.band_names)
lutf.load_lut(prod.lutfine, indband)
lutc.load_lut(prod.lutcoarse, indband)

# reproject lut array on the angles of the image
# angles are rounded to reduce the dims of interpolated LUT
sza_ = _utils.remove_na(np.unique(prod.sza.round(1)))
vza_ = _utils.remove_na(np.unique(prod.vza.round(1)))
azi_ = _utils.remove_na(np.unique(prod.razi.round(0)))
lutf.interp_n_slice(sza_, vza_, azi_)
lutc.interp_n_slice(sza_, vza_, azi_)
aotlut = np.array(lutf.aot, dtype=prod.type)

# plotting
lutf.refl.isel(sza=0, vza=0, azi=[0, 10]).plot(x='wl', hue='aot', col='azi')
lutc.refl.isel(sza=0, vza=0, azi=[0, 10]).plot(x='wl', hue='aot', col='azi')

##################################
# GET ANCILLARY DATA (AEROSOL)
##################################
aero = acutils.aerosol()
aot550rast = np.zeros([prod.height, prod.width], dtype=prod.type, order='F')
aotscarast = np.zeros([prod.height, prod.width], dtype=prod.type, order='F')
# ssarast = np.zeros([prod.N,prod.width, prod.height], dtype=prod.type)
aotrast = np.zeros([prod.N, prod.height, prod.width], dtype=prod.type, order='F')
fcoefrast = np.zeros([prod.height, prod.width], dtype=prod.type, order='F')
# set mean SSA value to adjust scattering AOT when info is not available
ssacoef = 0.99
if prod.aerosol == 'cds_forecast':

    cams_file = os.path.join(prod.cams_folder, prod.date.strftime('%Y'), prod.date.strftime('%Y-%m') +
                             '_month_cams-global-atmospheric-composition-forecasts.nc')
    prod.aux.get_xr_cams_cds_aerosol(cams_file, prod, lutf, lutc)
    aot550rast = prod.aux.aot_sca_550  # .T
    fcoefrast = prod.aux.fcoef
    aotrast = prod.aux.aot_grs
    aotscarast = prod.aux.aot_sca_grs
    logging.info(f'aot550rast shape {aot550rast.shape}')

else:
    # AERONET data
    if (prod.aerosol == 'aeronet'):
        prod.set_aeronetfile(aeronet_file)
        try:
            prod.aux.Aeronet.import_aeronet_data(aero, prod.aeronetfile, prod.date)
        except:
            logging.info('Error: No aeronet data in the +/- 2day time window.')
            sys.exit()

        prod.aux.aot_wl = aero.wavelengths
        prod.aux.aot = aero.aot
        aotscarast = ssacoef * prod.aux.aot
        prod.aux.aot550 = aero.aot550
        aot550rast.fill(prod.aux.aot550)

    # CAMS dataset
    elif (prod.aerosol == 'cams_forecast') | (prod.aerosol == 'cams_reanalysis'):

        # monthly file
        # target = Path(
        #     os.path.join(prod.cams_folder, prod.date.strftime('%Y'), prod.date.strftime('%Y-%m') + '_month_' +
        #                  prod.aerosol + '.nc'))
        cams_file = os.path.join(prod.cams_folder, prod.date.strftime('%Y'),
                                 prod.date.strftime('%Y-%m') + '_month_' +
                                 prod.aerosol + '.nc')
        prod.aux.get_xr_cams_aerosol(cams_file, prod.product)
        aotscarast = ssacoef * prod.aux.aot
        aot550rast = prod.aux.aot550rast  # .T
        logging.info(f'aot550rast shape {aot550rast.shape}')
    # CAMS new cds dataset (available from 26 June 2018 12UTC)

    elif (prod.aerosol == 'user_model'):
        prod.aux.aot550 = aot550
        prod.angstrom = angstrom
        prod.aux.aot = prod.aux.aot550 * (np.array(prod.wl) / 550) ** (-prod.angstrom)
        aotscarast = ssacoef * prod.aux.aot
        prod.aux.aot_wl = prod.wl
        aot550rast.fill(prod.aux.aot550)

    else:
        prod.aux.aot550 = 0.1
        prod.angstrom = 1
        prod.aux.aot = prod.aux.aot550 * (np.array(prod.wl) / 550) ** (-prod.angstrom)
        aotscarast = ssacoef * prod.aux.aot
        prod.aux.aot_wl = prod.wl
        aot550rast.fill(prod.aux.aot550)

        logging.info("No aerosol data provided, set to default: aot550=01, angstrom=1")

    # set spectral aot for satellite bands
    aero.fit_spectral_aot(prod.aux.aot_wl, prod.aux.aot)
    prod.aot = aero.get_spectral_aot(np.array(prod.wl))
    prod.aot550 = prod.aux.aot550

    # normalization of Cext to get spectral dependence of fine and coarse modes
    nCext_f = lutf.Cext / lutf.Cext550
    nCext_c = lutc.Cext / lutc.Cext550
    logging.info('param aerosol {nCext_f}, {nCext_c}, {prod.aot}')
    aero.fit_aero(nCext_f, nCext_c, prod.aot / prod.aot550)
    logging.info(f'{aero.fcoef} {aero.fcoef.astype(prod.type)}')
    prod.fcoef = aero.fcoef.astype(prod.type)
    fcoefrast.fill(prod.fcoef[0])
    for i in range(prod.N):
        aotrast[i].fill(prod.aot[i])

prod.rot = prod.sensordata.rot

####################################
#     Set SMAC parameters for
#    absorbing gases correction
####################################
logging.info('loading SMAC algorithm...')
smac = acutils.smac(prod.sensordata.smac_bands, prod.sensordata.smac_dir)
smac.set_gas_param()
smac.set_values(o3du=prod.aux.o3du, h2o=prod.aux.h2o)
smac.set_standard_values(prod.pressure_msl)
prod.aux.no2 = smac.uno2

######################################
# arrays allocation
# reshaping for fortran binding
######################################

aot550guess = np.zeros(prod.width, dtype=prod.type)
rtoaf = np.zeros((lutf.aot.__len__(), prod.N, prod.width), dtype=prod.type, order='F')
rtoac = np.zeros((lutc.aot.__len__(), prod.N, prod.width), dtype=prod.type, order='F')
maskpixels_ = np.full(prod.width, 1, dtype=prod.type, order='F')

w, h = prod.width, prod.height

rcorr = np.zeros((prod.N, h, w), dtype=prod.type)  # , order='F').T
rcorrg = np.zeros((prod.N, h, w), dtype=prod.type)  # , order='F').T
aot550pix = np.zeros((w, h), dtype=prod.type, order='F').T
betapix = np.zeros((w, h), dtype=prod.type, order='F').T
brdfpix = np.zeros((w, h), dtype=prod.type, order='F').T

prod.l2_product.getBand('SZA').writePixels(0, 0, w, h, prod.sza)
prod.l2_product.getBand('VZA').writePixels(0, 0, w, h, np.array(prod.vza[1]))
prod.l2_product.getBand('AZI').writePixels(0, 0, w, h, np.array(prod.razi[1]))

######################################
#      Add terrain attributes
######################################
if dem:
    sza_mean, sazi_mean = np.nanmean(prod.sza), np.nanmean(prod.sazi)
    prod.slope, prod.shade = _utils.get_dem_attributes(prod.elevation, sza=sza_mean, sun_azi=sazi_mean)
    # add elevation band
    prod.l2_product.getBand('elevation').writePixels(0, 0, w, h, prod.elevation)
    prod.l2_product.getBand('slope').writePixels(0, 0, w, h, prod.slope)
    prod.l2_product.getBand('shade').writePixels(0, 0, w, h, prod.shade)

######################################
#      MAIN LOOP
######################################
logging.info('processing ' + file + '...')

######################################
#      First step: AOT adjustment
######################################
# TODO slice/reshape.. to use the two SWIR bands only
# xblock, yblock = 25, 25
# TODO aot estimation on water pixels, check better spatial resolution

######################################
#      Second step: Atmosphere and surface correction
######################################
# TODO put chunck size in config yaml file
xblock, yblock = 512, 512
for iy in range(0, w, yblock):
    print('process row ' + str(iy) + ' / ' + str(h))
    yc = iy + yblock
    if yc > w:
        yc = w
    for ix in range(0, h, xblock):
        # print('process col ' + str(ix) + ' / ' + str(w))
        xc = ix + xblock
        if xc > h:
            xc = h
        sza = prod.sza[ix:xc, iy:yc]
        xshape, yshape = sza.shape

        if (xshape == 0) or (yshape == 0):
            continue

        razi = prod.razi[:, ix:xc, iy:yc]
        vza = prod.vza[:, ix:xc, iy:yc]
        muv = prod.muv[:, ix:xc, iy:yc]
        mu0 = prod.mu0[ix:xc, iy:yc]
        mask = prod.mask[ix:xc, iy:yc]
        flags = prod.flags[ix:xc, iy:yc]
        band_rad = prod.band_rad[:, ix:xc, iy:yc]

        maskpixels = maskpixels_
        if allpixels:
            maskpixels = maskpixels * 0
        elif waterdetect_only:
            # print('watermask', prod.watermask[ix:xc, iy:yc].shape)
            maskpixels[prod.watermask[ix:xc, iy:yc] == 1] = 0
        else:
            maskpixels = mask

        if dem:
            elev = prod.elevation[ix:xc, iy:yc]
            pressure = acutils.misc.get_pressure(elev, prod.pressure_msl)
            pressure_corr = pressure / prod.pressure_ref
        else:
            pressure_corr = np.full((xshape, yshape), prod.pressure / prod.pressure_ref)
        pressure_corr = np.array(pressure_corr, dtype=prod.type, order='F')

        # ---------
        # if maja L2A image provided, use AOT_MAJA product
        # AOT_maja seems to be largely overestimated, before further analyses: force usage of CAMS instead
        # if maja:
        #     aot550guess = prod.aot_maja[i]
        #     # aot550guess[aot550guess < 0.01] = 0.01
        # else:
        #     aot550guess = aot550rast[i]

        aot550guess = np.array(aot550rast[ix:xc, iy:yc], dtype=prod.type, order='F')
        fcoef = np.array(fcoefrast[ix:xc, iy:yc], dtype=prod.type, order='F')
        aot_tot = np.array(aotrast[:, ix:xc, iy:yc], dtype=prod.type, order='F')
        if prod.aerosol == 'cds_forecast':
            aot_sca = np.array(aotscarast[:, ix:xc, iy:yc], dtype=prod.type, order='F')
        else:
            aot_sca = aot_tot

        for iband in range(prod.N):
            # correct for gaseous absorption
            tg = smac.compute_gas_trans(iband, prod.pressure_msl, mu0.reshape(-1),
                                        muv[iband].reshape(-1)).reshape(xshape, yshape)
            band_rad[iband] = band_rad[iband] / tg

        p = grs_solver.grs.main_algo(xshape, yshape, *lutf.refl.shape,
                                     aotlut, sza_, azi_, vza_,
                                     lutf.refl, lutc.refl, lutf.Cext, lutc.Cext,
                                     vza, sza, razi, band_rad, maskpixels,
                                     prod.wl, pressure_corr, prod.sensordata.rg, prod.solar_irr, prod.rot,
                                     aot_tot, aot_sca, aot550guess, fcoef,
                                     prod.nodata, prod.rrs)

        rcorr[:, ix:xc, iy:yc] = p[0]  # np.transpose(p[0],(0,2,1))
        rcorrg[:, ix:xc, iy:yc] = p[1]  # np.transpose(p[1],(0,2,1))
        aot550pix[ix:xc, iy:yc] = p[2]
        brdfpix[ix:xc, iy:yc] = p[3]

        # TODO improve checksum scheme
        prod.checksum('row:col, ' + str(ix) + ':' + str(iy))

rcorr[rcorr == prod.nodata] = np.nan
rcorrg[rcorrg == prod.nodata] = np.nan

ndwi_corr = np.array((rcorrg[prod.sensordata.NDWI_vis] - rcorrg[prod.sensordata.NDWI_nir]) / \
                     (rcorrg[prod.sensordata.NDWI_vis] + rcorrg[prod.sensordata.NDWI_nir]))
# set flags
prod.flags = prod.flags + \
             ((prod.mask == 1) +
              (np.array((rcorr[1] < -0.01) | (rcorr[2] < -0.01)) << 1) +
              ((prod.mask == 2) << 2) +
              (((ndwi_corr < prod.sensordata.NDWI_threshold[0]) | (
                      ndwi_corr > prod.sensordata.NDWI_threshold[1])) << 3) +
              ((rcorrg[prod.sensordata.high_nir[0]] > prod.sensordata.high_nir[1]) << 4)
              )
print(w, h, prod.flags.astype(np.uint32).shape, brdfpix.shape, aot550pix.shape, rcorr.shape)
prod.l2_product.getBand('flags').writePixels(0, 0, w, h, np.array(prod.flags.astype(np.uint32)))
prod.l2_product.getBand('BRDFg').writePixels(0, 0, w, h, brdfpix)
prod.l2_product.getBand("aot550").writePixels(0, 0, w, h, aot550pix)
for iband in range(prod.N):
    prod.l2_product.getBand(prod.output + '_' + prod.band_names[iband]). \
        writePixels(0, 0, w, h, rcorr[iband])
    prod.l2_product.getBand(prod.output + '_g_' + prod.band_names[iband]). \
        writePixels(0, 0, w, h, rcorrg[iband])

prod.finalize_product()

if unzip:
    # remove unzipped files (Sentinel files)
    shutil.rmtree(file, ignore_errors=True)

if untar:
    # remove untared files (Landsat files)
    shutil.rmtree(tmp_dir, ignore_errors=True)
