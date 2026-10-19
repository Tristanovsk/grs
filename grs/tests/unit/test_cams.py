"""CAMS data (grs.CamsProduct) on the test subsets."""

import os

import pytest
import xarray as xr

from grs import CamsProduct


@pytest.fixture(scope='module')
def raster(l1c_path):
    with xr.open_dataset(l1c_path) as ds:
        return ds[['sza']].load()


def test_explicit_file(raster):
    cams = CamsProduct(raster, cams_file='/some/where/cams_forecast_2021-02.nc')
    assert (cams.dir, cams.file) == ('/some/where', 'cams_forecast_2021-02.nc')


def test_file_of_the_image(raster, tmp_path):
    cams = CamsProduct(raster, dir=str(tmp_path / 'cams'), suffix='_31TFJ')
    assert cams.filepath == os.path.join(str(tmp_path / 'cams'), '2021-02-11-forecast_31TFJ.nc')


def test_area(raster):
    cams = CamsProduct(raster, cams_file='cams.nc')
    north, west, south, east = cams.area                       # image footprint enlarged by 1 deg
    assert north > south and east > west
    assert 3 < west < 4 and 5 < east < 6


def test_load(raster, cams_path):
    cams = CamsProduct(raster, cams_file=cams_path)
    cams.load()
    for var in ('aod550', 'sp', 'gtco3', 'tcwv', 'u10', 'v10'):
        assert var in cams.raster, var
    assert float(cams.raster.aod550.mean()) > 0
