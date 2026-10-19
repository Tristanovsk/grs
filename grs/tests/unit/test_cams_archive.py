"""CAMS files: YYYY/MM/DD archive, file of an image, download into the archive."""

import os

import pytest
import xarray as xr

import grs.cams
from grs.cams import CAMS_DATASET, CamsProduct, cams_daily_file
from grs.grs_process import _cams_from_archive


@pytest.fixture(scope='module')
def raster(l1c_path):
    with xr.open_dataset(l1c_path) as ds:
        return ds[['sza']].load()


@pytest.fixture
def fake_ads(monkeypatch):
    """cdsapi client writing a small file instead of requesting the ADS; returns the requests."""
    calls = []

    class Client:
        def retrieve(self, name, request, target):
            calls.append(dict(name=name, request=request, target=target))
            with open(target, 'w') as f:
                f.write('area file')

    monkeypatch.setattr(grs.cams.cdsapi, 'Client', Client)
    return calls


def test_daily_file_path():
    path = cams_daily_file('2026-09-09', '/cams')
    assert path == f'/cams/2026/09/09/2026-09-09-{CAMS_DATASET}.nc'
    assert cams_daily_file('2026-09-09', '/cams', '_31TFJ').endswith(f'{CAMS_DATASET}_31TFJ.nc')


def test_file_of_the_area_when_no_global_file(raster, tmp_path):
    cams = CamsProduct(raster, dir=str(tmp_path), suffix='_31TFJ')
    assert cams.filepath == cams_daily_file('2021-02-11', str(tmp_path), '_31TFJ')
    assert not os.path.exists(cams.dir)               # nothing created before a download


def test_area_name_without_suffix(raster, tmp_path):
    cams = CamsProduct(raster, dir=str(tmp_path))
    assert '_N' in cams.file and '_E' in cams.file


def test_global_file_used_when_present(raster, tmp_path):
    global_file = cams_daily_file('2021-02-11', str(tmp_path))
    os.makedirs(os.path.dirname(global_file))
    open(global_file, 'w').close()
    assert CamsProduct(raster, dir=str(tmp_path), suffix='_31TFJ').filepath == global_file


def test_download_never_overwrites_global_file(raster, tmp_path, fake_ads):
    global_file = cams_daily_file('2021-02-11', str(tmp_path))
    os.makedirs(os.path.dirname(global_file))
    with open(global_file, 'w') as f:
        f.write('GLOBAL')
    cams = CamsProduct(raster, dir=str(tmp_path), suffix='_31TFJ')
    cams.cams_download()
    assert fake_ads[-1]['target'] == cams_daily_file('2021-02-11', str(tmp_path), '_31TFJ')
    assert fake_ads[-1]['request']['area'] == cams.area
    assert open(global_file).read() == 'GLOBAL'
    assert cams.filepath == fake_ads[-1]['target']


def test_former_layout_still_read(raster, tmp_path):
    raster = raster.assign_attrs(tile='31TFJ')
    legacy = tmp_path / '31TFJ' / '2021-02-11-forecast_31TFJ.nc'
    legacy.parent.mkdir()
    legacy.touch()
    assert _cams_from_archive(raster, str(tmp_path)).filepath == str(legacy)
    global_file = cams_daily_file('2021-02-11', str(tmp_path))
    os.makedirs(os.path.dirname(global_file))
    open(global_file, 'w').close()
    assert _cams_from_archive(raster, str(tmp_path)).filepath == global_file


def test_explicit_file(raster):
    cams = CamsProduct(raster, cams_file='/some/where/cams_forecast_2021-02.nc')
    assert (cams.dir, cams.file) == ('/some/where', 'cams_forecast_2021-02.nc')


def test_load_test_subset(raster, cams_path):
    cams = CamsProduct(raster, cams_file=cams_path)
    cams.load()
    for var in ('aod550', 'sp', 'gtco3', 'tcwv', 'wind'):
        assert var in cams.raster, var
    assert float(cams.raster.aod550.mean()) > 0
