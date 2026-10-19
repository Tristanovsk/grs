"""Input product (grs.Product) and aerosol parameterization (grs.Aerosol) on the test L1C subset."""

import numpy as np
import pytest
import xarray as xr

from grs import Aerosol, Product


@pytest.fixture(scope='module')
def prod(l1c_path):
    with xr.open_dataset(l1c_path) as ds:
        return Product(ds.load())


def test_product(prod):
    assert prod.sensor.startswith('S2')
    assert prod.raster.bands.dims == ('wl', 'y', 'x')
    assert (prod.width, prod.height) == (128, 128)
    assert np.issubdtype(prod.wl.dtype, np.integer)
    assert prod.raster.attrs['version'] == __import__('grs').__version__
    assert 0 < float(prod.sza_mean) < 90 and 0 <= float(prod.vza_mean) < 20
    assert 2 < float(prod.air_mass_mean) < 5


def test_product_bounds(prod):
    assert prod.lonmin < prod.lonmax and prod.latmin < prod.latmax
    assert 4 < prod.lonmin < 5 and 43 < prod.latmin < 44            # Camargue
    assert prod.xmax - prod.xmin == pytest.approx(128 * 60)


def test_spectral_aot_fit():
    wl = np.array([469., 550., 670., 865., 1240.])
    aot = 0.2 * (wl / 550) ** -1.3                                   # Angstrom law
    aero = Aerosol()
    aero.fit_spectral_aot(wl, aot)
    np.testing.assert_allclose(aero.get_spectral_aot(wl), aot, rtol=1e-6)
    assert aero.get_spectral_aot(np.array([550.]))[0] == pytest.approx(0.2, rel=1e-6)
