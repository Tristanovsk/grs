'''
Main program
'''

from pathlib import Path
from esasnappy import ProductData, ProductIO

import os, shutil
import zipfile
import tarfile
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import grs
from grs import config as cfg
from grs import acutils
from grs import AuxData

from grs.anglegen import *
from grs.fortran.grs import main_algo as grs_solver
from grs.fortran.grs_a import main_algo as grs_a_solver




file='/datalake/S2-L1C/21MYT/2022/01/14/S2A_MSIL1C_20220114T141051_N0301_R110_T21MYT_20220114T155152.SAFE'
outfile='./test.nc'
wkt=None; sensor=None; aerosol='default'; ancillary=None; altitude=0;
dem=True; aeronet_file=None; aot550=0.1; angstrom=1; resolution=None; unzip=False; untar=False;
startrow=0; allpixels=False; maja_xml=None; waterdetect_file=None; waterdetect_only=False;
memory_safe=False; angleonly=False; grs_a=False; output='Rrs'
'''
Main program calling all GRS steps

:param file: Input file to be processed
:param outfile: Absolute path of the output file
:param wkt: Well-Known-Text format defining the area of interest for which the image is subset
:param sensor: Set the sensor type: S2A, S2B, LANDSAT_5, LANDSAT_7, LANDSAT_8
            (by default sensor type is retrieved from input file name)
:param aerosol: aerosol data base to use within the processing
           DB: cams_forecast, cams_reanalysis, cds_forecast, aeronet, user_model, default
:param ancillary: if None, value is set to that of aerosol
:param altitude: provide altitude if `dem` is set as `False`
:param dem: if True digital elevation model is applied for per-pixel pressure calculation (data from SNAP/SRTM)
:param aeronet_file: optional aeronet file to be used for aerosol calculations
:param maja_xml: optional use of mask from MAJA L2A images, path to xml ID of the L2A image
:param waterdetect_file: optional use of water mask from waterdetect algorithm,
            path to the appropriate WaterDetect data file
:param waterdetect_only: if True and waterdetect file is provided, process only the pixels masked as "water"
:param resolution: pixel resolution in meters (integer)
:param unzip: if True input file is unzipped before processing,
              NB: unzipped files are removed at the end of the process
:param startrow: row number of the resampled and subset image on which the process starts, recommended value 0
                NB: this option is used to in the context of operational processing of massive dataset
:param allpixels: force to process all pixels even they are flagged as "Vegetation" or "Non-water"
:param angleonly: if true, grs is used to compute angle parameters only (no atmo correction is applied)
:param output: set the unit of the retrievals:

         * 'Lwn', normalized water-leaving radiance (in  :math:`mW cm^{-2} sr^{-1} \mu m^{-1})`

         * 'Rrs', remote sensing reflectance (in  :math:`sr^{-1}`)

         {default: 'Rrs'}
:param grs_a: switch to grs-a algorithm (Lwn and aerosol) if True

:return:
'''

##################################
# Get sensor auxiliary data
##################################

print('Get sensor auxiliary data')
_utils = grs.utils()
if sensor == None:
    sensor = _utils.get_sensor(file)
sensordata = AuxData.sensordata(sensor)
if resolution == None:
    resolution = sensordata.resolution
indband = sensordata.indband

if ancillary == None:
    ancillary = aerosol

##################################
# Read L1C product
##################################


file_orig = file
# unzip if needed
if unzip:
    print('unzipping...')
    tmpzip = zipfile.ZipFile(file)
    tmpzip.extractall(cfg.tmp_dir)
    file = os.path.join(cfg.tmp_dir, tmpzip.namelist()[0])
tartmp = None
if untar:
    basename = os.path.basename(file).replace('.tgz', '')
    basename = os.path.basename(basename).replace('.tar.gz', '')
    tmp_dir = os.path.join(cfg.tmp_dir, basename)
    # open tar archive to extract files for data loading
    tmpzip = tarfile.open(file)
    tmpzip.extractall(tmp_dir)
    file = glob.glob(os.path.join(tmp_dir, '*MTL.*'))[0]
    # open tar archive to add potential file (e.g., angle files) - InvalidHeaderError
    if not any(['solar' in f for f in glob.glob(os.path.join(tmp_dir, '*'))]):
        tartmp = tarfile.open(os.path.join(cfg.tmp_dir, os.path.basename(file_orig)), 'w:gz')

print("Reading...")
print(file)
product = ProductIO.readProduct(file)

##################################
# Generate l2h object
##################################
l2h = grs.info(product, sensordata, aerosol, ancillary, output)
l2h.headerfile = file

##################################
# GET METADATA
##################################
print('getting metadata...')
# TODO clean up this part and other metadata to be loaded
if 'S2' in sensor:
    meta = l2h.Product.getMetadataRoot().getElement('Level-1C_User_Product').getElement(
        'General_Info').getElement(
        'Product_Image_Characteristics').getElement('Reflectance_Conversion')
    l2h.U = float(str(meta.getAttribute('U').getData()))
    l2h.solar_irr = np.zeros(len(indband), dtype=np.float32)
    for i, iband in zip(range(len(indband)), indband):
        l2h.solar_irr[i] = float(str(meta.getElement('Solar_Irradiance_List').getAttributeAt(iband).getData()))

else:
    meta = l2h.Product.getMetadataRoot().getElement("L1_METADATA_FILE").getElement("IMAGE_ATTRIBUTES")
    l2h.U = float(str(meta.getAttribute('EARTH_SUN_DISTANCE').getData())) ** 2
    l2h.solar_irr = np.array(l2h.SensorData.solar_irr)[indband]

# convert into mW cm-2 um-1
l2h.solar_irr = l2h.solar_irr / 10

##################################
# GENERATE BAND ANGLES (LANDSAT)
##################################
anggen = False
if 'LANDSAT_8' in sensor:
    anggen = angle_generator().landsat(l2h)
elif 'LANDSAT' in sensor:
    anggen = angle_generator().landsat_tm(l2h)
#print('anggen = {}'.format(anggen))
if anggen:
    if tartmp:
        print('writing angles to input file: ' + file_orig)
        # TODO finalize this part to add angle files to original tar.gz image (e.g., LC8*.tgz)
        # copy tgz image and add angle files
        shutil.move(file_orig, os.path.join(os.path.dirname(file_orig), 'saves', os.path.basename(file_orig)))
        for filename in os.listdir(tmp_dir):
            tartmp.add(os.path.join(tmp_dir, filename), filename)
        tartmp.close()
        shutil.move(os.path.join(cfg.tmp_dir, os.path.basename(file_orig)), file_orig)


##################################
# RESAMPLE TO A UNIQUE RESOLUTION
##################################
print('resampling...')
if 'S2' in sensor:
    if memory_safe:
        l2h.Product = _utils.generic_resampler(l2h.Product, resolution=resolution)  # , method='Nearest')
    else:
        l2h.Product = _utils.s2_resampler(l2h.Product, resolution=resolution)
else:
    l2h.Product = _utils.resampler(l2h.Product, resolution=resolution)  # , upmethod='Nearest')

##################################
# SUBSET TO AREA OF INTEREST
##################################
print('subsetting...')
try:
    if wkt is not None:
        l2h.Product = _utils.get_subset(l2h.Product, wkt)
    l2h.get_product_info()
except:
    if unzip:
        # remove unzipped files (Sentinel files)
        shutil.rmtree(file, ignore_errors=True)
    if untar:
        # remove untared files (Landsat files)
        shutil.rmtree(tmp_dir, ignore_errors=True)
    raise NameError('No data available for requested area')

l2h.set_outfile(outfile)
l2h.wkt, lonmin, lonmax, latmin, latmax = _utils.get_extent(l2h.Product)
l2h.crs = str(l2h.Product.getBand(l2h.band_names[0]).getGeoCoding().getImageCRS())

##################################
# Fetch optional mask products
# resample for common resolution
# subset to ROI
##################################
print('fetching falgs...')
maja, waterdetect = None, None
if maja_xml:
    try:
        maja = ProductIO.readProduct(maja_xml)
        maja = _utils.resampler(maja, resolution=resolution)
        maja = _utils.get_subset(maja, wkt)

    except:
        print('!!! issues with ' + maja_xml + '; please check if file exists')
        raise

if waterdetect_file:
    try:
        waterdetect = ProductIO.readProduct(waterdetect_file)
        waterdetect = _utils.get_subset(waterdetect, wkt)
        waterdetect = _utils.resampler(waterdetect, resolution=resolution)

    except:
        print('!!! issues with ' + waterdetect_file + '; please check if file exists')
        raise

##################################
## ADD ELEVATION BAND
##################################
print('adding elevation band...')
if dem:
    print('add elevation band')
    high_latitude = (latmax >= 60) | (latmin <= -60)
    l2h.get_elevation(high_latitude)

##################################
# GET IMAGE AND RASTER PROPERTIES
##################################
print('load raster data...')
l2h.get_bands(l2h.band_names[0:2])
l2h.print_info()

##################################
# SET NEW BAND FOR MASKING
##################################
print('add ndwi mask')
l2h.Product.addBand('ndwi', ProductData.TYPE_FLOAT32)
l2h.ndwi_band = l2h.Product.getBand('ndwi')
l2h.ndwi_band.ensureRasterData()
l2h.ndwi_band.loadRasterData()

##################################
# GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
##################################
print('getting CAMS data...')
l2h.aux = AuxData.cams()

if ancillary != 'default':
    if l2h.Aerosol == 'cds_forecast':
        target = os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'), l2h.date.strftime('%Y-%m') +
                              '_month_cams-global-atmospheric-composition-forecasts.nc')
        l2h.aux.get_cams_ancillary(target, l2h.date, l2h.wkt, param=['msl', 'gtco3', 'tcwv', 'tcno2', 't2m'])
    else:
        target = Path(os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'),
                                   l2h.date.strftime('%Y-%m') + '_month_' + l2h.ancillary + '.nc'))
        # do not load here since already implemented elsewhere in CNES HPC
        # l2h.aux.load_cams_data(target, l2h.date, data_type=l2h.ancillary)
        l2h.aux.get_cams_ancillary(target, l2h.date, l2h.wkt)

## uncomment this part to use ecmwf files provided in the .SAFE format
# if 'S2' in sensor:
#     l2h.aux.get_tile_dir(file)
#     l2h.aux.get_aux_dir()
#     l2h.aux.get_ecmwf_data()

# get pressure at the scene altitude
l2h.pressure_msl = l2h.aux.msl  # acutils.misc.get_pressure(altitude, l2h.aux.msl)
if dem:
    altitude = l2h.elevation
    altitude[altitude < -200] = 0
l2h.pressure = acutils.Misc.get_pressure(altitude, l2h.pressure_msl)  # l2h.aux.pressure = l2h.pressure

######################################
#      Create output l2 product
#          'l2_product'
######################################
print('creating L2 output product')
l2h.create_product(maja=maja, waterdetect=waterdetect)
try:
    l2h.load_data()
except:
    if unzip:
        # remove unzipped files (Sentinel files)
        shutil.rmtree(file, ignore_errors=True)
    if untar:
        # remove untared files (Landsat files)
        shutil.rmtree(tmp_dir, ignore_errors=True)
    raise NameError('No data available for requested area')
l2h.load_flags()

#####################################
# LOAD LUT FOR ATMOSPHERIC CORRECTION
#####################################
print('loading lut...', l2h.lutfine)
lutf = acutils.lut(l2h.band_names)
lutc = acutils.lut(l2h.band_names)
lutf.load_lut(l2h.lutfine, indband)
lutc.load_lut(l2h.lutcoarse, indband)
# reproject lut array on the angles of the image
# angles are rounded to reduce the dims of interpolated LUT
sza_ = _utils.remove_na(np.unique(l2h.sza.round(1)))
vza_ = _utils.remove_na(np.unique(l2h.vza.round(1)))
azi_ = _utils.remove_na(np.unique(l2h.razi.round(0)))
lutf.interp_n_slice(sza_,vza_,azi_)
lutc.interp_n_slice(sza_,vza_,azi_)
aotlut = np.array(lutf.aot, dtype=l2h._type)

##################################
# GET ANCILLARY DATA (AEROSOL)
##################################
aero = acutils.Aerosol()
aot550rast = np.zeros([l2h.height, l2h.width], dtype=l2h._type, order='F')
aotscarast = np.zeros([l2h.height, l2h.width], dtype=l2h._type, order='F')
# ssarast = np.zeros([l2h.N,l2h.width, l2h.height], dtype=l2h.type)
aotrast = np.zeros([l2h.N, l2h.height, l2h.width], dtype=l2h._type, order='F')
fcoefrast = np.zeros([l2h.height, l2h.width], dtype=l2h._type, order='F')
# set mean SSA value to adjust scattering AOT when info is not available
ssacoef = 0.99
if l2h.Aerosol == 'cds_forecast':

    cams_file = os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'), l2h.date.strftime('%Y-%m') +
                             '_month_cams-global-atmospheric-composition-forecasts.nc')
    l2h.aux.get_xr_cams_cds_aerosol(cams_file, l2h, lutf, lutc)
    aot550rast = l2h.aux.aot_sca_550  # .T
    fcoefrast = l2h.aux.fcoef
    aotrast = l2h.aux.aot_grs
    aotscarast = l2h.aux.aot_sca_grs
    print('aot550rast shape', aot550rast.shape)

else:
    # AERONET data
    if (l2h.Aerosol == 'aeronet'):
        l2h.set_aeronetfile(aeronet_file)
        try:
            l2h.aux.Aeronet.import_aeronet_data(aero, l2h.aeronetfile, l2h.date)
        except:
            print('Error: No aeronet data in the +/- 2day time window.')
            sys.exit()

        l2h.aux.aot_wl = aero.wavelengths
        l2h.aux.aot = aero.aot
        aotscarast = ssacoef*l2h.aux.aot
        l2h.aux.aot550 = aero.aot550
        aot550rast.fill(l2h.aux.aot550)

    # CAMS dataset
    elif (l2h.Aerosol == 'cams_forecast') | (l2h.Aerosol == 'cams_reanalysis'):

        # monthly file
        # target = Path(
        #     os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'), l2h.date.strftime('%Y-%m') + '_month_' +
        #                  l2h.aerosol + '.nc'))
        cams_file = os.path.join(l2h.cams_folder, l2h.date.strftime('%Y'),
                                 l2h.date.strftime('%Y-%m') + '_month_' +
                                 l2h.Aerosol + '.nc')
        l2h.aux.get_xr_cams_aerosol(cams_file, l2h.Product)
        aotscarast = ssacoef*l2h.aux.aot
        aot550rast = l2h.aux.aot550rast  # .T
        print('aot550rast shape', aot550rast.shape)
    # CAMS new cds dataset (available from 26 June 2018 12UTC)

    elif (l2h.Aerosol == 'user_model'):
        l2h.aux.aot550 = aot550
        l2h.angstrom = angstrom
        l2h.aux.aot = l2h.aux.aot550 * (np.array(l2h.wl) / 550) ** (-l2h.angstrom)
        aotscarast = ssacoef*l2h.aux.aot
        l2h.aux.aot_wl = l2h.wl
        aot550rast.fill(l2h.aux.aot550)

    else:
        l2h.aux.aot550 = 0.1
        l2h.angstrom = 1
        l2h.aux.aot = l2h.aux.aot550 * (np.array(l2h.wl) / 550) ** (-l2h.angstrom)
        aotscarast = ssacoef*l2h.aux.aot
        l2h.aux.aot_wl = l2h.wl
        aot550rast.fill(l2h.aux.aot550)

        print("No aerosol data provided, set to default: aot550=01, angstrom=1")

    # set spectral aot for satellite bands
    aero.fit_spectral_aot(l2h.aux.aot_wl, l2h.aux.aot)
    l2h.aot = aero.get_spectral_aot(np.array(l2h.wl))
    l2h.aot550 = l2h.aux.aot550

    # normalization of Cext to get spectral dependence of fine and coarse modes
    nCext_f = lutf.Cext / lutf.Cext550
    nCext_c = lutc.Cext / lutc.Cext550
    print('param aerosol', nCext_f, nCext_c, l2h.aot)
    aero.fit_aero(nCext_f, nCext_c, l2h.aot / l2h.aot550)
    print(aero.fcoef, aero.fcoef.astype(l2h._type))
    l2h.fcoef = aero.fcoef.astype(l2h._type)
    fcoefrast.fill(l2h.fcoef[0])
    for i in range(l2h.N):
        aotrast[i].fill(l2h.aot[i])

l2h.rot = l2h.SensorData.rot

####################################
#     Set SMAC parameters for
#    absorbing gases correction
####################################
print('loading SMAC algorithm...')
smac = acutils.smac(l2h.SensorData.smac_bands, l2h.SensorData.smac_dir)
smac.set_gas_param()
smac.set_values(o3du=l2h.aux.o3du, h2o=l2h.aux.h2o)
smac.set_standard_values(l2h.pressure_msl)
l2h.aux.no2 = smac.uno2



######################################
# arrays allocation
# reshaping for fortran binding
######################################

aotlut = np.array(lutf.aot, dtype=l2h._type, order='F')
vzalut = np.array(lutf.vza, dtype=l2h._type, order='F')
szalut = np.array(lutf.sza, dtype=l2h._type, order='F')
razilut = np.array(lutf.azi, dtype=l2h._type, order='F')
rlut_f = np.array(lutf.refl, dtype=l2h._type, order='F')
rlut_c = np.array(lutc.refl, dtype=l2h._type, order='F')
grid_lut = (szalut, razilut, vzalut)

aot550guess = np.zeros(l2h.width, dtype=l2h._type)
rtoaf = np.zeros((lutf.aot.__len__(), l2h.N, l2h.width), dtype=l2h._type, order='F')
rtoac = np.zeros((lutc.aot.__len__(), l2h.N, l2h.width), dtype=l2h._type, order='F')
maskpixels_ = np.full(l2h.width, 1, dtype=l2h._type, order='F')

w, h = l2h.width, l2h.height

l2h.l2_product.getBand('SZA').writePixels(0, 0, w, h, l2h.sza)
l2h.l2_product.getBand('VZA').writePixels(0, 0, w, h, np.array(l2h.vza[1]))
l2h.l2_product.getBand('AZI').writePixels(0, 0, w, h, np.array(l2h.razi[1]))

if dem:
    # add elevation band
    l2h.l2_product.getBand('elevation').writePixels(0, 0, w, h, l2h.elevation)

######################################
#      MAIN LOOP
######################################
print('processing '+file+'...')
for i in range(startrow, l2h.height):
    # print('process row ' + str(i) + ' / ' + str(l2h.height))

    sza = l2h.sza[i]
    razi = l2h.razi[:, i]
    vza = l2h.vza[:, i]
    muv = l2h.muv[:, i]
    mu0 = l2h.mu0[i]
    mask = l2h.mask[i]
    flags = l2h.flags[i]
    band_rad = l2h.band_rad[:, i]

    maskpixels = maskpixels_
    if allpixels:
        maskpixels = maskpixels * 0
    elif waterdetect_only:
        maskpixels[l2h.watermask[i] == 1] = 0
    else:
        maskpixels = mask

    if dem:
        elev = l2h.elevation[i]
        pressure = acutils.Misc.get_pressure(elev, l2h.pressure_msl)
        pressure_corr = pressure / l2h.pressure_ref
    else:
        pressure_corr = [l2h.pressure / l2h.pressure_ref] * l2h.width
    pressure_corr = np.array(pressure_corr, dtype=l2h._type, order='F')

    # ---------
    # if maja L2A image provided, use AOT_MAJA product
    # AOT_maja seems to be largely overestimated, before further analyses: force usage of CAMS instead
    # if maja:
    #     aot550guess = l2h.aot_maja[i]
    #     # aot550guess[aot550guess < 0.01] = 0.01
    # else:
    #     aot550guess = aot550rast[i]

    aot550guess = np.array(aot550rast[i])
    fcoef = np.array(fcoefrast[i])
    aot_tot = np.array(aotrast[:, i])
    if l2h.Aerosol == 'cds_forecast':
        aot_sca = np.array(aotscarast[:, i])
    else:
        aot_sca = aot_tot

    for iband in range(l2h.N):
        # preparing lut data
        grid_pix = list(zip(sza, razi[iband], vza[iband]))

        for iaot in range(aotlut.__len__()):
            rtoaf[iaot, iband] = lutf.interp_lut(grid_lut, rlut_f[iband][iaot, ...], grid_pix)
            rtoac[iaot, iband] = lutc.interp_lut(grid_lut, rlut_c[iband][iaot, ...], grid_pix)

        # correct for gaseous absorption
        tg = smac.compute_gas_trans(iband, l2h.pressure_msl, mu0, muv[iband])
        band_rad[iband] = band_rad[iband] / tg

    if grs_a:
        rcorr, rcorrg, aot550pix, brdfpix = grs_a_solver.main_algo(l2h.width, l2h.N, aotlut.__len__(),
                                                                   vza, sza, razi, band_rad, maskpixels, l2h.wl,
                                                                   pressure_corr, aotlut, rtoaf, rtoac,
                                                                   lutf.Cext, lutc.Cext, lutf.Cext550,
                                                                   lutc.Cext550,
                                                                   l2h.SensorData.rg, l2h.solar_irr, l2h.rot,
                                                                   aot_tot,
                                                                   aot550guess, fcoef, l2h.nodata, l2h.rrs)
    else:

        rcorr, rcorrg, aot550pix, brdfpix = grs_solver.main_algo(l2h.width, l2h.N, aotlut.__len__(),
                                                                 vza, sza, razi, band_rad, maskpixels, l2h.wl,
                                                                 pressure_corr, aotlut, rtoaf, rtoac, lutf.Cext,
                                                                 lutc.Cext,
                                                                 l2h.SensorData.rg, l2h.solar_irr, l2h.rot,
                                                                 aot_tot, aot_sca, aot550guess, fcoef, l2h.nodata,
                                                                 l2h.rrs)

    # reshape for snap modules
    rcorr[rcorr == l2h.nodata] = np.nan
    rcorrg[rcorrg == l2h.nodata] = np.nan
    # rcorr = np.ma.array(rcorr.T, mask=rcorr.T == l2h.nodata, fill_value=np.nan)  # .tolist()
    # rcorrg = np.ma.array(rcorrg.T, mask=rcorrg.T == l2h.nodata, fill_value=np.nan)  # .tolist()

    ndwi_corr = np.array((rcorrg[l2h.SensorData.NDWI_vis] - rcorrg[l2h.SensorData.NDWI_nir]) / \
                         (rcorrg[l2h.SensorData.NDWI_vis] + rcorrg[l2h.SensorData.NDWI_nir]))
    # set flags

    flags = flags + \
            ((mask == 1) +
             (np.array((rcorr[1] < -0.01) | (rcorr[2] < -0.01)) << 1) +
             ((mask == 2) << 2) +
             (((ndwi_corr < l2h.SensorData.NDWI_threshold[0]) | (
                     ndwi_corr > l2h.SensorData.NDWI_threshold[1])) << 3) +
             ((rcorrg[l2h.SensorData.high_nir[0]] > l2h.SensorData.high_nir[1]) << 4)
             )

    for iband in range(l2h.N):
        l2h.l2_product.getBand(l2h.output + '_' + l2h.band_names[iband]). \
            writePixels(0, i, l2h.width, 1, np.array(rcorr[iband]))
        l2h.l2_product.getBand(l2h.output + '_g_' + l2h.band_names[iband]). \
            writePixels(0, i, l2h.width, 1, np.array(rcorrg[iband]))

    l2h.l2_product.getBand('flags').writePixels(0, i, l2h.width, 1, flags.astype(np.uint32))
    l2h.l2_product.getBand('BRDFg').writePixels(0, i, l2h.width, 1, brdfpix)
    l2h.l2_product.getBand("aot550").writePixels(0, i, l2h.width, 1, aot550pix)

    # TODO improve checksum scheme
    l2h.checksum('row ' + str(i)+ ' / ' + str(l2h.height))

l2h.finalize_product()

if unzip:
    # remove unzipped files (Sentinel files)
    shutil.rmtree(file, ignore_errors=True)

if untar:
    # remove untared files (Landsat files)
    shutil.rmtree(tmp_dir, ignore_errors=True)
