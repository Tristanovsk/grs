"""Early stop of Process.execute for the images without water (min_water_area)."""

import types

import numpy as np
import xarray as xr

from grs.grs_process import Process, _water_area


def _prod(water_pixels, size=50, pixel=20.):
    '''Product with `water_pixels` water pixels (blue > SWIR1, low SWIR2) on a land image.'''
    wl = [490, 1610, 2190]
    bands = np.empty((3, size, size), dtype=np.float32)
    bands[:] = np.array([0.05, 0.2, 0.1], dtype=np.float32)[:, None, None]    # land: SWIR1 > blue
    flat = bands.reshape(3, -1)
    flat[:, :water_pixels] = np.array([0.05, 0.01, 0.005], dtype=np.float32)[:, None]   # water
    bands[0, -1, -1] = np.nan    # nodata pixel: never water
    raster = xr.Dataset({'bands': (('wl', 'y', 'x'), bands)},
                        coords={'wl': wl, 'x': np.arange(size) * pixel, 'y': -np.arange(size) * pixel})
    return types.SimpleNamespace(raster=raster, bvis=490, bswir=1610, bswir2=2190,
                                 vis_swir_index_threshold=0., sunglint_threshold=0.2)


def test_water_area():
    # 20 m pixels: 400 m² each
    assert _water_area(_prod(0)) == 0
    assert np.isclose(_water_area(_prod(100)), 100 * 400 / 1e6)


def test_write_output_without_product():
    process = Process()
    process.no_water = True
    assert process.write_output() is None
