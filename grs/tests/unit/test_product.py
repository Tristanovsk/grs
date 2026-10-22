"""Level-1 product (grs.product.Product) of the test L1C subset."""

import numpy as np
import pytest
import xarray as xr

from grs.product import Product


@pytest.fixture(scope='module')
def prod(l1c_path):
    return Product(xr.open_dataset(l1c_path))


def test_metadata(prod):
    assert prod.sensor == 'S2A'
    assert prod.date_str == '2021-02-11T10:41:51.024Z'
    assert prod.date.year == 2021 and prod.date.month == 2 and prod.date.day == 11


def test_grid(prod):
    assert (prod.width, prod.height) == (128, 128)
    # 128 pixels at 60 m in UTM 31N (tile 31TFJ)
    assert (prod.xmin, prod.ymin, prod.xmax, prod.ymax) == pytest.approx((618540, 4812300, 626220, 4819980))
    assert prod.xmax - prod.xmin == pytest.approx(128 * 60)
    assert prod.raster.rio.crs.to_epsg() == 32631


def test_geographic_bounds(prod):
    # the Camargue, around the Vaccarès lagoon
    assert prod.lonmin == pytest.approx(4.4652, abs=1e-3)
    assert prod.lonmax == pytest.approx(4.5618, abs=1e-3)
    assert prod.latmin == pytest.approx(43.4530, abs=1e-3)
    assert prod.latmax == pytest.approx(43.5234, abs=1e-3)


def test_spectral_response(prod):
    srf = prod.raster.SRF
    assert srf.sizes['wl'] == 13                            # the 13 bands of MSI
    np.testing.assert_allclose(srf.max('wl_hr'), 1, rtol=1e-6)  # normalized to 1 at the peak
    assert np.nanmin(srf.values) >= 0


def test_gas_luts_packaged(prod):
    # gas absorption LUTs shipped with the package (not downloaded)
    assert 'wl' in prod.gas_lut.coords
    assert 'Twv' in prod.Twv_lut
    assert prod.air_mass_mean == pytest.approx(2.971, abs=1e-3)
