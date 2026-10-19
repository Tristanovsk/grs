"""Gaseous absorption correction, in memory and block by block (dask)."""

import numpy as np
import pytest
import xarray as xr

from grs.grs_process import _correct_gas_absorption

WL = [443., 490.]


def _bands():
    data = np.random.default_rng(0).random((2, 100, 100)).astype(np.float32)
    return xr.DataArray(data, dims=('wl', 'y', 'x'),
                        coords=dict(wl=WL, y=np.arange(100.)[::-1], x=np.arange(100.), spatial_ref=0))


@pytest.mark.parametrize('spatial_ref', [0, 32631])
def test_dask_equals_numpy(spatial_ref):
    # transmittance on a coarse grid with its own spatial_ref (as the CAMS raster): the dask path failed
    # when it differed from the one of the bands
    tg = xr.DataArray(np.random.default_rng(1).random((2, 4, 4)) * 0.1 + 0.9, dims=('wl', 'y', 'x'),
                      coords=dict(wl=WL, y=np.linspace(110, -10, 4), x=np.linspace(-10, 110, 4),
                                  spatial_ref=spatial_ref))
    lazy = _correct_gas_absorption(_bands().chunk({'y': 50, 'x': 50, 'wl': -1}), tg).compute()
    eager = _correct_gas_absorption(_bands(), tg)
    np.testing.assert_allclose(lazy.values, eager.values, rtol=1e-6)
    assert int(lazy.spatial_ref) == 0
