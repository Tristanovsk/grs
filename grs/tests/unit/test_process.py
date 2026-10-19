"""End-to-end test: L1C subset -> L2A NetCDF product (needs the look-up tables, see conftest.py)."""

import numpy as np
import pytest
import xarray as xr

import grs

EXPECTED_VARIABLES = ('Rrs', 'BRDFg', 'aot550', 'flags', 'mask', 'sza', 'vza', 'raa')


@pytest.fixture(scope='module')
def l2a(l2a_files):
    with xr.open_dataset(l2a_files[0]) as ds:
        yield ds.load()


@pytest.fixture(scope='module')
def l1c(l1c_path):
    with xr.open_dataset(l1c_path) as ds:
        yield ds.load()


def test_variables_and_grid(l2a, l1c):
    missing = [v for v in EXPECTED_VARIABLES if v not in l2a]
    assert not missing
    np.testing.assert_allclose(l2a.x, l1c.x)
    np.testing.assert_allclose(l2a.y, l1c.y)
    assert l2a.Rrs.dims == ('wl', 'y', 'x')
    assert l2a.sizes['wl'] == 11


def test_attributes(l2a):
    assert l2a.attrs['version'] == grs.__version__
    assert l2a.attrs['tile'] == '31TFJ'


def test_ancillary_file(l2a_files):
    with xr.open_dataset(l2a_files[1]) as anc:
        for var in ('aod550', 'sp', 'gtco3', 'tcwv'):
            assert var in anc, var


def test_mask(l2a):
    assert set(np.unique(l2a.mask.values)) <= {0, 1}
    assert 0.2 < float((l2a.mask == 0).mean()) < 0.95      # water, land and shoreline in the subset


def test_rrs_on_water(l2a):
    water = (l2a.mask == 0).values
    for wl in (443, 560, 665, 865):
        rrs = l2a.Rrs.sel(wl=wl).values
        assert np.isfinite(rrs[water]).mean() > 0.99, wl
    assert 0 < np.median(l2a.Rrs.sel(wl=560).values[water]) < 0.1


def test_aot(l2a):
    aot = l2a.aot550.values
    assert np.all(np.isfinite(aot))
    assert np.all(aot >= 0) and float(aot.max()) < 1.5
