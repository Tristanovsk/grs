"""Gaseous transmittance (grs.utils.GaseousTransmittance) and pressure (grs.utils.Misc)."""

import numpy as np
import pytest
import xarray as xr

from grs import CamsProduct, utils
from grs.product import Product


@pytest.fixture(scope='module')
def gas(l1c_path, cams_path):
    prod = Product(xr.open_dataset(l1c_path))
    cams = CamsProduct(prod.raster, cams_file=cams_path)
    cams.load()
    return utils.GaseousTransmittance(prod, cams)


def test_inputs(gas):
    assert (gas.xmin, gas.ymin, gas.xmax, gas.ymax) == pytest.approx((618540, 4812300, 626220, 4819980))
    # surface pressure of CAMS in hPa
    assert float(gas.pressure.mean()) == pytest.approx(1014.2, abs=1)


def test_total_transmittance(gas):
    tg = gas.get_gaseous_transmittance()
    assert dict(tg.sizes) == {'y': 12, 'x': 12, 'wl': 13}
    assert np.all((tg.values > 0) & (tg.values <= 1))
    by_band = tg.mean(['x', 'y'])
    # water vapour bands strongly absorbed: B10 (1375 nm, the most) and B9 (945 nm);
    # B8A (865 nm) almost transparent
    assert float(by_band.idxmin('wl')) == pytest.approx(1375, abs=10)
    assert float(by_band.sel(wl=1375, method='nearest')) < 0.05
    assert float(by_band.sel(wl=945, method='nearest')) < 0.6
    assert float(by_band.sel(wl=865, method='nearest')) > 0.99


def test_background_transmittance(gas):
    tb = gas.Tgas_background()
    assert dict(tb.sizes) == {'y': 12, 'x': 12, 'wl': 13}
    assert np.all((tb.values > 0.97) & (tb.values <= 1))


@pytest.mark.parametrize('alt, psl, expected', [
    (0., 1013.25, 1013.25),
    (1000., 998., 885.236756238),
])
def test_pressure_at_altitude(alt, psl, expected):
    assert utils.Misc.get_pressure(alt, psl) == pytest.approx(expected, abs=1e-8)
