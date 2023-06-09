import glob
import os
import xml.etree.ElementTree as ET

import numpy as np
from numba import jit

import pandas as pd
import geopandas as gpd
import xarray as xr
import xmltodict

from rasterio.features import rasterize
import scipy.odr as odr
from affine import Affine
from osgeo import gdal, ogr
import cartopy.crs as ccrs
from pyproj import CRS
import eoreader as eo
from eoreader.reader import Reader

opj = os.path.join
BAND_NAMES = np.array(['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'])
BAND_NAMES_EOREADER = np.array(['CA', 'BLUE', 'GREEN', 'RED', 'VRE_1',
                                'VRE_2', 'VRE_3', 'NARROW_NIR', 'NIR',
                                'WV', 'SWIR_CIRRUS', 'SWIR_1', 'SWIR_2'])

BAND_ID = [b.replace('B', '') for b in BAND_NAMES]
NATIVE_RESOLUTION = [60, 10, 10, 10, 20, 20, 20, 10, 20, 60, 60, 20, 20]
WAVELENGTH = np.array([443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190])
BAND_WIDTH = [20, 65, 35, 30, 15, 15, 20, 115, 20, 20, 30, 90, 180]

INFO = pd.DataFrame({'bandId': range(13),
                     'ESA': BAND_NAMES,
                     'EOREADER': BAND_NAMES_EOREADER,
                     'Wavelength (nm)': WAVELENGTH,
                     'Band width (nm)': BAND_WIDTH,
                     'Resolution (m)': NATIVE_RESOLUTION}).set_index('bandId').T


class s2image():
    def __init__(self, imageSAFE, band_idx=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                 resolution=20, verbose=False, **kwargs):

        abspath = os.path.abspath(imageSAFE)
        dirroot, basename = os.path.split(abspath)
        self.verbose = verbose
        self.band_idx = band_idx
        self.resolution = resolution
        self.INFO = INFO[band_idx]

        # --------------------------------
        # define interpolation parameters
        # --------------------------------
        # tile for 10m resolution: width,height = 10980,10980
        # tile for 20m resolution: width,height = 5490,5490
        # tile for 60m resolution: width,height = 1830,1830
        if resolution == 10:
            self.width, self.height = 10980, 10980
        elif resolution == 20:
            self.width, self.height = 5490, 5490
        elif resolution == 60:
            self.width, self.height = 1830, 1830

        self.xml_granule = glob.glob(opj(imageSAFE, 'GRANULE', '*', 'MTD_TL.xml'))[0]
        self.xml_file = glob.glob(opj(imageSAFE, 'MTD*.xml'))[0]

        # get metadata and angle file
        ds = gdal.Open(self.xml_file)
        self.metadata = ds.GetMetadata()
        self.metadata2 = xmltodict.parse(ds.GetMetadata('xml:SENTINEL2')[0])
        __ = []
        for _ in self.metadata2['n1:Level-1C_User_Product']['n1:General_Info']['Product_Image_Characteristics'][
            'Reflectance_Conversion']['Solar_Irradiance_List']['SOLAR_IRRADIANCE']:
            __.append([int(_['@bandId']), float(_['#text'])])
        self.solar_irradiance = np.array(__)

        # Spectral Response Functions
        SRFs = []
        wl_hr = np.arange(400, 2350)
        for _ in self.metadata2['n1:Level-1C_User_Product'][
            'n1:General_Info']['Product_Image_Characteristics']['Spectral_Information_List']['Spectral_Information']:

            bandID = int(_['@bandId'])
            if not self.band_idx.__contains__(bandID):
                continue

            wl_min, wl_max = float(_['Wavelength']['MIN']['#text']), float(_['Wavelength']['MAX']['#text'])
            step = float(_['Spectral_Response']['STEP']['#text'])
            wl = np.arange(wl_min, wl_max + step, step)
            RSF = np.asarray(_['Spectral_Response']['VALUES'].split(), dtype=np.float32)
            SRFs.append(xr.DataArray(RSF, coords=dict(wl_hr=wl), name='SRF').interp(wl_hr=wl_hr).assign_coords(
                dict(wl=WAVELENGTH[bandID])))
        self.SRFs = xr.concat(SRFs, dim='wl')
        self.SRFs.attrs['description'] = 'Spectral response function of each band'

        # Open instance of eoreader
        reader = Reader()

        # Open the product
        prod = reader.open(imageSAFE, remove_tmp=True, **kwargs)
        self.prod = prod
        self.processing_baseline = prod._processing_baseline
        self.datetime = prod.datetime

        # save geographic data
        self.extent = prod.extent()
        self.bounds = self.extent.bounds
        minx, miny, maxx, maxy = self.bounds.values[0]
        self.crs = self.prod.crs()
        self.epsg = self.extent.crs.to_epsg()
        str_epsg = str(self.epsg)
        zone = str_epsg[-2:]
        is_south = str_epsg[2] == 7
        self.proj = ccrs.UTM(zone, is_south)
        # self.proj = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': is_south})
        self.transform = Affine(resolution, 0., minx, 0., -resolution, maxy)

        # -------------------------
        # interpolation
        # -------------------------
        # set the indexing of numpy and xarray
        self.indexing = 'xy'
        self.new_x = np.linspace(minx, maxx, self.width)
        self.new_y = np.linspace(miny, maxy, self.height)[::-1]

        # ---------------------------------
        # getting appropriate version of mask driver
        # ---------------------------------
        # TODO generalize to load all available masks
        if self.processing_baseline < 4:
            self._open_mask = prod._open_mask_lt_4_0
        else:
            self._open_mask = prod._open_mask_gt_4_0

    def load_product(self, add_time=False, **kwargs):

        self.load_bands(add_time=False, **kwargs)
        self.load_geom()
        self.prod.attrs = self.metadata
        self.prod.attrs['satellite'] = self.prod.attrs['PRODUCT_URI'].split('_')[0]
        self.prod.attrs['solar_irradiance'] = self.solar_irradiance[:, 1]
        self.prod.attrs['solar_irradiance_unit'] = 'W/m²/µm'
        self.prod.attrs['acquisition_date'] = self.prod.attrs['DATATAKE_1_DATATAKE_SENSING_START']
        # add spectral response function
        self.prod = xr.merge([self.prod, self.SRFs]).drop_vars('bandID')

    def load_bands(self, add_time=False, **kwargs):

        # ----------------------------------
        # getting bands
        # ----------------------------------
        bands = self.prod.stack(list(BAND_NAMES_EOREADER[self.band_idx]), resolution=self.resolution, **kwargs)
        bands = bands.rename({'z': 'bandID'})

        # ----------------------------------
        # setting up coordinates and dimensions
        # ----------------------------------
        self.prod = bands.assign_coords(wl=('bandID', self.INFO.loc['Wavelength (nm)'])). \
            swap_dims({'bandID': 'wl'}).drop({'band', 'bandID', 'variable'})
        self.prod = self.prod.assign_coords(bandID=('wl', self.INFO.loc['ESA'].values))
        self.prod = self.prod.to_dataset(name='bands', promote_attrs=True)

        # add time
        if add_time:
            self.prod = self.prod.assign_coords(time=self.datetime).expand_dims('time')
        # self.prod.clear()

    @staticmethod
    def parse_angular_grid_node(node):
        '''
        Internal parsing function for angular grids
        :return:
        '''
        values = []
        for c in node.find('Values_List'):
            values.append(np.array([float(t) for t in c.text.split()]))
            values_array = np.stack(values)
        return values_array

    @staticmethod
    def set_crs(arr, crs):
        arr.rio.set_crs(crs, inplace=True)
        arr.rio.write_crs(inplace=True)

    def get_raw_angles(self):

        minx, miny, maxx, maxy = self.bounds.values[0]
        with open(self.xml_granule) as xml_file:
            tree = ET.parse(xml_file)
            root = tree.getroot()

        raw_sza = self.parse_angular_grid_node(root.find('.//Tile_Angles/Sun_Angles_Grid/Zenith'))
        raw_sazi = self.parse_angular_grid_node(root.find('.//Tile_Angles/Sun_Angles_Grid/Azimuth'))

        # compute x and y for angle grids
        # check dimension
        Nx, Ny = raw_sza.shape

        xang = np.linspace(minx, maxx, Nx)
        yang = np.linspace(miny, maxy, Ny)[::-1]

        raw_sun_ang = xr.Dataset(data_vars=dict(sza=(['y', 'x'], raw_sza),
                                                sazi=(['y', 'x'], raw_sazi)),
                                 coords=dict(x=xang, y=yang))
        self.set_crs(raw_sun_ang, self.crs)
        self.raw_sun_ang = raw_sun_ang

        # ---------------------------------
        # getting viewing geometry datacube
        # ---------------------------------
        bandIds, detectorIds = [], []
        for angleID in root.findall('.//Tile_Angles/Viewing_Incidence_Angles_Grids'):
            bandIds.append(int(angleID.attrib['bandId']))
            detectorIds.append(int(angleID.attrib['detectorId']))
        Nband, Ndetector = np.max(bandIds) + 1, np.max(detectorIds) + 1

        # allocate/fill rasters
        vza, vazi = np.full((Nband, Ndetector, Nx, Ny), np.nan, dtype=float), np.full((Nband, Ndetector, Nx, Ny),
                                                                                      np.nan,
                                                                                      dtype=float)

        for angleID in root.findall('.//Tile_Angles/Viewing_Incidence_Angles_Grids'):
            iband = int(angleID.attrib['bandId'])
            idetector = int(angleID.attrib['detectorId'])
            vza[iband, idetector] = self.parse_angular_grid_node(angleID.find('Zenith'))
            vazi[iband, idetector] = self.parse_angular_grid_node(angleID.find('Azimuth'))

        raw_view_ang = xr.Dataset(data_vars=dict(vza=(['bandId', 'detectorId', 'y', 'x'], vza),
                                                 vazi=(['bandId', 'detectorId', 'y', 'x'], vazi)),
                                  coords=dict(bandId=range(Nband),
                                              detectorId=range(Ndetector),
                                              x=xang, y=yang))
        self.set_crs(raw_view_ang, self.crs)

        # clean up Dataset (remove empty slices)
        raw_view_ang = raw_view_ang.dropna('detectorId', how='all')
        self.raw_view_ang = raw_view_ang

        # set number of different detectors
        self.detector_num = len(raw_view_ang.detectorId)

        return

    def load_geom(self, method='linear'):

        self.get_raw_angles()
        self.get_all_band_angles(method=method)

    @staticmethod
    def linfit(beta, x):
        return beta[0] * x[0] + beta[1] * x[1] + beta[2]

    @staticmethod
    @jit(nopython=True)  # "uint16[:,:](float64[:],float64[:],float64[:,:],float64[:,:],intp,intp)",
    def lin2D(arr, x, y, mask, betas, detector_offset=0, scale_factor=100):

        Nx, Ny = mask.shape

        for ii in range(Nx):
            for jj in range(Ny):
                detect = mask[ii, jj]
                if detect == 0:
                    continue
                beta = betas[detect - detector_offset]
                val = beta[0] * x[jj] + beta[1] * y[ii] + beta[2]
                # compression using simple int8 and scale factor
                arr[ii, jj] = (val * scale_factor)

    def data_fitting(self, x0, y0, arr):
        # ---------------------------------
        # ODR multilinear regression
        # ---------------------------------
        xgrid, ygrid = np.meshgrid(x0, y0, indexing=self.indexing)

        # vectorize
        values = arr.values.flatten()
        x_ = xgrid.flatten()
        y_ = ygrid.flatten()

        # remove NaN
        idx = ~np.isnan(values)
        values = values[idx]
        points = np.empty((2, len(values)))
        points[0] = x_[idx]
        points[1] = y_[idx]

        # set ODR fitting
        mean = np.nanmean(values)
        linear = odr.Model(self.linfit)
        data = odr.Data(points, values)
        beta0 = [0, 0, mean]

        # proceed with ODR fitting
        fit = odr.ODR(data, linear, beta0=beta0)
        resfit = fit.run()

        if self.verbose:
            resfit.pprint()

        return resfit.beta

    def get_detector_mask(self, bandId=0, resolution=20, detector_mask_name='DETFOO'):

        if self.processing_baseline < 4:
            mask_df = self._open_mask(detector_mask_name, BAND_ID[bandId])
            detector_num = mask_df.gml_id.str.split('-', expand=True).values[:, 2]
            poly_shp = [[geom, int(value)] for geom, value in zip(mask_df.geometry, detector_num)]

            mask = rasterize(shapes=poly_shp,
                             out_shape=(self.height, self.width),
                             transform=self.transform)
        else:

            mask = self._open_mask(detector_mask_name, BAND_ID[bandId], resolution=resolution).astype(np.int8)
            mask = mask.squeeze()
        return np.array(mask)

    def get_band_angle_as_numpy(self, xarr, bandId=0, resolution=20,
                                detector_mask_name='DETFOO', compress=False,
                                ):
        '''

        :param xarr:
        :param bandId:
        :param resolution:
        :param detector_mask_name:

        :return:
        '''

        detector_offset = xarr.detectorId.values.min()
        mask = self.get_detector_mask(bandId=bandId, resolution=resolution)

        # TODO check how to avoid taking the nodata value "0" when coarsening the raster
        # TODO for the moment this induces bad detector number at the edge of the image swath
        # mask = _open_mask(detector_mask_name, BAND_ID[bandId], resolution=NATIVE_RESOLUTION[bandId])
        # # mask nodata value
        # mask = mask.where(mask!=0)
        # # resample mask at the desired resolution
        # if resolution != NATIVE_RESOLUTION[iband]:
        #     mask = mask.interp(x=new_x, y=new_y, method='nearest')
        # # compress mask into int8
        # mask = mask.astype(np.int8)

        x, y = self.new_x, self.new_y
        # self.prod.clear()
        betas = np.full((self.detector_num, 3), np.nan)
        xarr_ = xarr.sel(bandId=bandId)
        for id in range(self.detector_num):
            # --------------------------------------------------------------
            # Linear 2D-fitting to get the function of the regression plane
            # --------------------------------------------------------------
            arr = xarr_.isel(detectorId=id).dropna('y', how='all').dropna('x', how='all')
            x0, y0 = arr.x.values, arr.y.values
            betas[id, :] = self.data_fitting(x0, y0, arr)

        if compress:
            # TODO experimental, needs to be tested
            # compression in uint16 (NB: range 0-65535)
            new_arr = np.full((self.width, self.height), np.nan, dtype=np.int16)
            # compute angles from betas values for each detector and band
            self.lin2D(new_arr, x, y, mask, betas, detector_offset=detector_offset, scale_factor=100)
        else:
            # compression in uint16 (NB: range 0-65535)
            new_arr = np.full((self.width, self.height), np.nan, dtype=np.float32)
            # compute angles from betas values for each detector and band
            self.lin2D(new_arr, x, y, mask, betas, detector_offset=detector_offset, scale_factor=1)

        del mask
        return new_arr

    @staticmethod
    def scat_angle(sza, vza, azi):
        '''
        self.azi: azimuth in rad for convention azi=180 when sun-sensor in opposition
        :return: scattering angle in deg
        '''

        sza = np.radians(sza)
        vza = np.radians(vza)
        azi = np.radians(azi)
        ang = -np.cos(sza) * np.cos(vza) - np.sin(sza) * np.sin(vza) * np.cos(azi)
        ang = np.arccos(ang)
        return np.degrees(ang)

    def get_all_band_angles(self, method='linear'):

        new_x, new_y = self.new_x, self.new_y
        band_idx = self.band_idx
        # -----------------------------------------------------------------
        # Sun angles (easy!) based on standard bidimensional interpolation
        # -----------------------------------------------------------------
        new_sun_ang = self.raw_sun_ang.interp(x=new_x, y=new_y, method=method)

        # ------------------------------------------------------
        # Viewing angles (not easy!) based on 2D-plane fitting
        # ------------------------------------------------------
        raw_vza = self.raw_view_ang.vza
        raw_vazi = self.raw_view_ang.vazi

        # ---------------------------------------------------------
        # convert vza, azi angles into cartesian vector components
        # (NOT NEEDED FOR THE MOMENT!!)
        # ---------------------------------------------------------
        # np.tan(np.deg2rad(view_ang.vza)) * np.sin(np.deg2rad(view_ang.vazi))
        # np.tan(np.deg2rad(view_ang.vza)) * np.cos(np.deg2rad(view_ang.vazi))
        # dx.name = 'dx'
        # dy.name = 'dy'

        new_vza, new_vazi = [], []
        for ibandId, bandId in enumerate(band_idx):
            if self.verbose:
                print('Band number ' + str(bandId) + ' is being loaded')
            new_vza.append(self.get_band_angle_as_numpy(raw_vza, bandId=bandId, resolution=self.resolution))
            new_vazi.append(self.get_band_angle_as_numpy(raw_vazi, bandId=bandId, resolution=self.resolution))
        raa = (np.array(new_vazi) - new_sun_ang.sazi.values) % 360

        self.prod['vza'] = xr.DataArray(np.array(new_vza), dims=['wl', 'y', 'x'])
        self.prod['raa'] = xr.DataArray(raa, dims=['wl', 'y', 'x'])
        self.prod['sza'] = xr.DataArray(new_sun_ang.sza.values, dims=['y', 'x'])
        #
        # xr.Dataset(data_vars=dict(vza=(['wl', 'y', 'x'], )),
        #                      coords=dict(wl=self.INFO.loc['Wavelength (nm)'],
        #
        # vazi = xr.Dataset(data_vars=dict(vazi=(['wl', 'y', 'x'], np.array(new_vazi))),
        #                   coords=dict(wl=self.INFO.loc['Wavelength (nm)'], x=new_x, y=new_y))
        # new_ang = xr.Dataset(data_vars=dict(vza=(['wl', 'y', 'x'], np.array(new_vza))),
        #                      coords=dict(wl=self.INFO.loc['Wavelength (nm)'], x=new_x, y=new_y))
        #
        # new_ang['sza'] = new_sun_ang.sza
        # new_ang['razi'] = vazi.vazi - new_sun_ang.sazi) % 360
        # # new_ang['scat_ang'] = scat_angle(new_ang.sza, new_ang.vza, new_ang.razi)

        # self.set_crs(new_ang, self.crs)

        del new_vza, new_vazi, raa, new_sun_ang

        return

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
