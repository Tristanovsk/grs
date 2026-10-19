"""Fixtures of the unit tests: small input data of `grs/tests/data` (see make_test_data.py)."""

import os

import pytest

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
# 128 x 128 pixels at 60 m of S2A 2021-02-11 over the Camargue (tile 31TFJ), and its CAMS data
L1C_SUBSET = os.path.join(DATA_DIR, 'S2A_MSIL1C_20210211_T31TFJ_camargue_subset.nc')
CAMS_SUBSET = os.path.join(DATA_DIR, 'cams_forecast_2021-02-11_camargue.nc')
# full Level-1 image of the subset, used only when it is on the computer
SAFE = ('/data/satellite/Sentinel-2/L1C/31TFJ/2021/02/11/'
        'S2A_MSIL1C_20210211T104151_N0500_R008_T31TFJ_20230608T124831.SAFE')


@pytest.fixture(scope='session')
def l1c_path():
    return L1C_SUBSET


@pytest.fixture(scope='session')
def cams_path():
    return CAMS_SUBSET


@pytest.fixture
def isolated_config(tmp_path, monkeypatch):
    """Configuration without the user file of the computer (empty user file, cache cleared)."""
    from grs.config import load_config
    user_file = tmp_path / 'config.yml'
    user_file.write_text('')
    monkeypatch.setenv('GRS_CONFIG', str(user_file))
    monkeypatch.delenv('GRS_NO_DOWNLOAD', raising=False)
    load_config.cache_clear()
    yield user_file
    load_config.cache_clear()


def luts_available():
    from grs.config import lut_status
    return all(exists for _, _, exists in lut_status())


requires_luts = pytest.mark.skipif(not luts_available(),
                                   reason='look-up tables not installed (grs_config --download)')


@pytest.fixture(scope='session')
def l2a_zarr(tmp_path_factory, l1c_path, cams_path):
    """L2A Zarr product made once by the full processing chain from the test L1C subset."""
    if not luts_available():
        pytest.skip('look-up tables not installed (grs_config --download)')
    import grs
    odir = str(tmp_path_factory.mktemp('l2a') / 'l2a.zarr')
    process = grs.Process()
    process.execute(l1c_path, cams_file=cams_path, resolution=60, odir=odir)
    process.write_output()
    return odir
