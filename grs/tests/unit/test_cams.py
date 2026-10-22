"""CAMS data (grs.cams.CamsProduct) on the grid of the test L1C subset."""

import numpy as np
import pytest
import xarray as xr

from grs import CamsProduct
from grs.product import Product


@pytest.fixture(scope='module')
def prod(l1c_path):
    return Product(xr.open_dataset(l1c_path))


@pytest.fixture(scope='module')
def cams(prod, cams_path):
    cams = CamsProduct(prod.raster, cams_file=cams_path)
    cams.load()
    return cams


def test_variables(cams):
    expected = {'u10', 'v10', 't2m', 'msl', 'sp', 'tcwv', 'gtco3', 'tcno2', 'tcco', 'tc_ch4',
                'aod469', 'aod550', 'aod670', 'aod865', 'aod1240'}
    assert expected <= set(cams.variables)
    assert expected <= set(cams.raster.data_vars)


def test_grid_covers_the_image(cams, prod):
    # coarse grid of the CAMS data spanning the image (same corners as the image grid)
    assert cams.raster.x.values[[0, -1]] == pytest.approx([prod.xmin, prod.xmax])
    assert cams.raster.y.values[[0, -1]] == pytest.approx([prod.ymax, prod.ymin])
    assert cams.raster.rio.crs == prod.raster.rio.crs


def test_values_of_the_day(cams):
    # winter morning over the Camargue: mistral (northerly wind) and about 1014 hPa
    assert float(cams.raster.v10.mean()) == pytest.approx(-6.2, abs=0.5)
    assert float(cams.raster.t2m.mean()) == pytest.approx(281.86, abs=0.1)
    assert float(cams.raster.sp.mean()) == pytest.approx(101420, abs=100)
    assert float(cams.raster.msl.mean()) == pytest.approx(101707, abs=20)
    assert float(cams.raster.tcwv.mean()) == pytest.approx(10.5, abs=0.5)


def test_aod_spectrum(cams):
    aod = cams.cams_aod
    assert aod.sizes['wl'] == 5
    assert np.all(aod.values > 0)
    # aerosol optical depth decreasing with the wavelength
    spectrum = aod.mean(['x', 'y']).sortby('wl').values
    assert np.all(np.diff(spectrum) < 0)
