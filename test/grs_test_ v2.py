'''
Main program
'''

import os

import matplotlib

matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr
import logging

from grs.drivers import driver_S2_SAFE as S2
from grs import Product, acutils, CamsProduct, L2aProduct
from grs.fortran.grs import main_algo as grs_solver
from grs import __version__
print(__version__)
opj = os.path.join
# from grs import product, l2a, utils, acutils, auxdata

# from grs.fortran.grs import main_algo as grs_solver
odir = '/data/satellite/sentinel2/L2A'
file = '/data/satellite/sentinel2/L1C/31TFJ/S2A_MSIL1C_20201004T104031_N0209_R008_T31TFJ_20201004T125253.SAFE'
file = '/media/harmel/vol1/Dropbox/satellite/S2/L1C/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE'
cams_file = '/media/harmel/vol1/Dropbox/satellite/S2/cnes/CAMS/2022-07-31-cams-global-atmospheric-composition-forecasts.nc'

file = '/data/satellite/S2/L1C//S2A_MSIL1C_20220713T103041_N0400_R108_T31TFJ_20220713T141110.SAFE'
cams_file = '/data/satellite/S2/cnes/CAMS/2022-07-13-cams-global-atmospheric-composition-forecasts.nc'
file_nc = file.replace('.SAFE', '.nc')
basename = os.path.basename(file)
ofile = basename.replace('.SAFE', '.nc').replace('L1C', 'L2Agrs')

bandIds = range(13)
resolution = 60

if not os.path.exists(file_nc):
    l1c = S2.s2image(file, band_idx=bandIds, resolution=resolution)
    print(l1c.crs)
    l1c.load_product()
    prod = Product(l1c)
    encoding = {
        'bands': {'dtype': 'int16', 'scale_factor': 0.00001, 'add_offset': .3, '_FillValue': -32768, "zlib": True,
                  "complevel": 6},
        'vza': {'dtype': 'int16', 'scale_factor': 0.01, '_FillValue': -32768, "zlib": True, "complevel": 6},
        'raa': {'dtype': 'int16', 'scale_factor': 0.01, 'add_offset': 180, '_FillValue': -32768, "zlib": True,
                "complevel": 6},
        'sza': {'dtype': 'int16', 'scale_factor': 0.01, 'add_offset': 30, '_FillValue': -32768, "zlib": True,
                "complevel": 6}}
    l1c.prod.to_netcdf(file_nc, encoding=encoding)
    l1c.prod.close()

else:
    prod = Product(xr.open_dataset(file_nc))

##################################
# Fetch optional mask products
# resample for common resolution
# subset to ROI
##################################
# TODO (if necessary) activate MAJA reader and flag retrieval
# prod.get_flags()

##################################
# GET ANCILLARY DATA (Pressure, O3, water vapor, NO2...
##################################
# prod.get_cams()
cams_dir = '/media/harmel/vol1/Dropbox/satellite/S2/cnes/CAMS'
cams = CamsProduct(prod, cams_file=cams_file)
cams.plot_params()

##################################
## ADD ELEVATION AND PRESSURE BAND
##################################
prod.get_elevation()
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


####################################
#    absorbing gases correction
####################################
gas_trans = acutils.GaseousTransmittance(prod, cams)
Tg_raster = gas_trans.get_gaseous_transmittance()
Tg_raster_coarse = gas_trans.Tg_tot_coarse

prod.raster['bands'] = prod.raster.bands / Tg_raster
prod.raster.bands.attrs['gas_absorption_correction'] = True

plt.figure()
Tg_raster.mean('x').mean('y').squeeze().plot(x='wl')
# Tg_raster.isel(wl=1).plot()
p = Tg_raster.plot.imshow(col='wl', col_wrap=3, robust=True, cmap=plt.cm.Spectral_r,
                          subplot_kws=dict(xlabel='', ylabel='', xticks=[], yticks=[]))

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

prod.raster['ndwi'] = ndwi
prod.raster.ndwi.attrs = {
    'description': 'Normalized difference spectral index between bands at ' + str(prod.b565) + ' and ' + str(
        prod.b865) + ' nm', 'units': '-'}
prod.raster['ndwi_swir'] = ndwi_swir
prod.raster.ndwi_swir.attrs = {
    'description': 'Normalized difference spectral index between bands at ' + str(prod.b565) + ' and ' + str(
        prod.b1600) + ' nm', 'units': '-'}

self = prod
masked_raster = prod.raster.bands.where(ndwi > self.ndwi_threshold). \
    where(b2200 < self.sunglint_threshold). \
    where(ndwi_swir > self.green_swir_index_threshold)

wl_process = [443, 490, 560, 665, 705,
              740, 783, 842, 865, 1610, 2190]
raster = masked_raster.sel(wl=wl_process)

Nx = prod.width
Ny = prod.height

aotlut = lutf.aot


def rounding(xarr, resol=1):
    vals = np.unique(xarr.round(resol))
    return vals[~np.isnan(vals)]


sza_ = rounding(prod.raster.sza, 1)
azi_ = rounding((180 - prod.raster.raa) % 360, 0)
vza_ = rounding(prod.raster.vza, 1)

aotlut = np.array(lutf.aot, dtype=prod._type)

# load OPAC LUT
lut_file = '/DATA/git/satellite_app/hgrs/data/lut/opac_osoaa_lut_v3.nc'
aero_lut = xr.open_dataset(lut_file)
aero_lut['wl']=aero_lut['wl']*1000

models = [ 'COAV_rh70', 'DESE_rh70']#'ANTA' + rh, 'ARCT' + rh,
wl = prod.central_wavelength
lut_aod = aero_lut.aot.sel(model=models, aot_ref=1).interp(wl=wl)
idx = np.abs((cams_aod / cams.aod550) - lut_aod).sum('wl').argmin()
opac_model = aero_lut.sel(model=models).model.values[idx]

lut_ = lutf.refl.sel(vza=slice(np.min(vza_)-1,np.max(vza_)+1),azi=slice(np.min(azi_)-1,np.max(azi_)+1))
mu0=np.linspace(0.01,1)
plt.figure()
for wl in range(8):
    z=np.polyfit(np.cos(np.radians(lutf.sza)),lutf.refl.isel(vza=2,azi=2,aot=2,wl=wl),2)
    p=np.poly1d(z)
    plt.plot(np.cos(np.radians(lutf.sza)),lutf.refl.isel(vza=2,azi=2,aot=2,wl=wl))
    plt.plot(mu0,p(mu0),ls='--',color='gray')

fine_refl = lutf.refl.interp(vza=vza_).interp(azi=azi_).interp(sza=sza_)
coarse_refl = lutc.refl.interp(vza=vza_).interp(azi=azi_).interp(sza=sza_)
lut_shape = fine_refl.shape
fine_Cext = lutf.Cext
coarse_Cext = lutc.Cext
vza = prod.raster.sel(wl=wl_process).vza.values
sza = prod.raster.sel(wl=wl_process).sza.values
razi = prod.raster.sel(wl=wl_process).raa.values
band_rad = raster.values
maskpixels = np.zeros((prod.height, prod.width))
wl_sat = wl_process
pressure_corr = cams.raster.sp.interp(x=raster.x, y=raster.y) * 1e-2 / prod.pressure_ref
eps_sunglint = prod.sensordata.rg
solar_irr = prod.solar_irradiance.sel(wl=wl_process).values
rot = prod.sensordata.rot

aot_tot_cams_res = cams.cams_aod.interp(wavelength=wl_process)
aot_sca_cams_res = aot_tot_cams_res * cams.cams_ssa.interp(wavelength=wl_process)

aot_tot = aot_tot_cams_res.interp(x=raster.x, y=raster.y)
aot_sca = aot_sca_cams_res.interp(x=raster.x, y=raster.y)
aot550guess = cams.raster.aod550.interp(x=raster.x, y=raster.y)
fcoef = np.full((prod.height, prod.width), 0.5)

rrs = prod.rrs

p = grs_solver.grs.main_algo(Nx, Ny, *lut_shape,
                             aotlut, sza_, azi_, vza_,
                             fine_refl, coarse_refl, fine_Cext, coarse_Cext,
                             vza, sza, razi, band_rad, maskpixels,
                             wl_sat, pressure_corr, eps_sunglint, solar_irr, rot,
                             aot_tot, aot_sca, aot550guess, fcoef, rrs)

rcorr, rcorrg, aot550pix, brdfpix = p

Rrs = xr.DataArray(rcorr, coords=raster.coords, name='Rrs')
Rrs_g = xr.DataArray(rcorrg, coords=raster.coords, name='Rrs_g')
aot550 = xr.DataArray(aot550pix, coords={'y': raster.y, 'x': raster.x}, name='aot550')
brdfg = xr.DataArray(brdfpix, coords={'y': raster.y, 'x': raster.x}, name='BRDFg')
l2_prod = xr.merge([Rrs, Rrs_g, aot550, brdfg])

l2a = L2aProduct(prod, l2_prod, cams, gas_trans)
l2a.to_netcdf(opj(odir, ofile))

##########################################
# PLOTTING SECTION
##########################################
plt.figure()
l2a.l2_prod.BRDFg.plot.imshow(robust=True, cmap=plt.cm.gray)

l2a.l2_prod.Rrs.plot.imshow(col='wl', col_wrap=4, robust=True, vmin=0, vmax=0.015, cmap=plt.cm.Spectral_r)

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
