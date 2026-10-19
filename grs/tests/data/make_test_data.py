"""Make the small input data of the unit tests of `grs/tests/unit` (stored in this folder).

- **L1C subset** (NetCDF): 128 x 128 pixels at 60 m (7.7 km) of the Sentinel-2A image of 2021-02-11
  over the Camargue lagoons and the coast (tile 31TFJ), read with GRSdriver: clear to very turbid
  waters, land and shoreline, as the test product of GRSl2bgen.
- **CAMS subset** (NetCDF): forecasts of the day of the image over the area of the image (+/- 1.5 deg).

Usage::

    python grs/tests/data/make_test_data.py [<L1C .SAFE>] [<monthly CAMS file>]
"""

import os
import sys
import warnings

import numpy as np
import xarray as xr

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
SAFE = ('/data/satellite/Sentinel-2/L1C/31TFJ/2021/02/11/'
        'S2A_MSIL1C_20210211T104151_N0500_R008_T31TFJ_20230608T124831.SAFE')
CAMS = '/data/cams/world/cams_forecast_2021-02.nc'
L1C_SUBSET = os.path.join(HERE, 'S2A_MSIL1C_20210211_T31TFJ_camargue_subset.nc')
CAMS_SUBSET = os.path.join(HERE, 'cams_forecast_2021-02-11_camargue.nc')
RESOLUTION = 60
SIZE = 128                                   # pixels
# lower left corner on the 60 m grid of the tile (upper left corner 600000, 4900020), bounds inside by 1 m so
# that only the SIZE x SIZE pixels of the window are read
X0, Y0 = 600000 + 309 * RESOLUTION, 4900020 - 1462 * RESOLUTION
BOUNDS = (X0 + 1, Y0 + 1, X0 + SIZE * RESOLUTION - 1, Y0 + SIZE * RESOLUTION - 1)   # EPSG:32631


def _clean_attrs(ds):
    '''Attributes that NetCDF can store (booleans as integers, None dropped).'''
    def clean(attrs):
        out = {}
        for k, v in attrs.items():
            if v is None:
                continue
            out[k] = int(v) if isinstance(v, (bool, np.bool_)) else v
        return out
    ds.attrs = clean(ds.attrs)
    for var in ds.variables.values():
        var.attrs = clean(var.attrs)
    return ds


def make_l1c(safe=SAFE):
    import GRSdriver
    l1c = GRSdriver.Sentinel2Driver(safe, resolution=RESOLUTION, subset=[BOUNDS])
    l1c.load_product()
    prod = l1c.prod.load()
    assert prod.sizes['x'] == SIZE and prod.sizes['y'] == SIZE, dict(prod.sizes)
    prod = _clean_attrs(prod)
    prod.attrs['test_subset'] = f'{os.path.basename(safe)}, bounds {BOUNDS} (EPSG:32631), {RESOLUTION} m'
    encoding = {v: {'zlib': True, 'complevel': 6} for v in prod.data_vars
                if prod[v].dtype.kind in 'fiu' and prod[v].ndim >= 2}
    for v in ('bands', 'vza', 'raa', 'sza'):     # float32 is enough for the tests
        if v in prod:
            prod[v] = prod[v].astype('float32')
    prod.to_netcdf(L1C_SUBSET, encoding=encoding)
    return prod


def make_cams(prod, cams=CAMS):
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(32631, 4326, always_xy=True).transform(
        [BOUNDS[0], BOUNDS[2]], [BOUNDS[1], BOUNDS[3]])
    ds = xr.open_dataset(cams)
    day = str(prod.time.dt.strftime('%Y-%m-%d').values)
    if 'time' in ds.dims:
        ds = ds.sel(time=day)
    ds = ds.sel(latitude=slice(max(lat) + 1.5, min(lat) - 1.5), longitude=slice(min(lon) - 1.5, max(lon) + 1.5))
    _clean_attrs(ds).load().to_netcdf(CAMS_SUBSET)
    return ds


if __name__ == '__main__':
    prod = make_l1c(*sys.argv[1:2])
    make_cams(prod, *sys.argv[2:3])
    for f in (L1C_SUBSET, CAMS_SUBSET):
        print(f'written {f} ({os.path.getsize(f) / 1e6:.2f} MB)')
