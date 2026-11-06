"""Fixtures of the unit tests: small input data of `grs/tests/data` (see make_test_data.py)."""

import os

import pytest

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
# 128 x 128 pixels at 60 m of S2A 2021-02-11 over the Camargue (tile 31TFJ), and its CAMS data
L1C_SUBSET = os.path.join(DATA_DIR, 'S2A_MSIL1C_20210211_T31TFJ_camargue_subset.nc')
CAMS_SUBSET = os.path.join(DATA_DIR, 'cams_forecast_2021-02-11_camargue.nc')


@pytest.fixture(scope='session')
def l1c_path():
    return L1C_SUBSET


@pytest.fixture(scope='session')
def cams_path():
    return CAMS_SUBSET


def lut_files():
    """
    Look-up tables used by the end-to-end tests: those of grs/config.yml, in the folder given by the
    environment variable GRS_TEST_GRSDATA if set (e.g. a local copy of grsdata).
    """
    from grs.grs_process import GRSDATA, TOALUT, TRANSLUT
    grsdata = os.environ.get('GRS_TEST_GRSDATA', GRSDATA)
    return os.path.join(grsdata, TOALUT), os.path.join(grsdata, TRANSLUT)


def luts_available():
    return all(os.path.exists(f) for f in lut_files())


@pytest.fixture(scope='session')
def l2a_files(tmp_path_factory, l1c_path, cams_path):
    """
    L2A product made once by the full processing chain from the test L1C subset: NetCDF files
    <odir>/<name>.nc and <odir>/<name>_anc.nc (ancillary data).
    """
    if not luts_available():
        pytest.skip('look-up tables not found (set GRS_TEST_GRSDATA to the folder of the LUTs)')
    import grs
    odir = tmp_path_factory.mktemp('l2a') / 'l2a'
    process = grs.Process()
    process.lut_file, process.trans_lut_file = lut_files()
    process.execute(l1c_path, cams_file=cams_path, resolution=60, odir=str(odir))
    process.write_output()
    return str(odir / 'l2a.nc'), str(odir / 'l2a_anc.nc')
