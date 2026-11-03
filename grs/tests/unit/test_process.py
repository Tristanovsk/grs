"""End-to-end test: L1C subset -> L2A Zarr product (needs the look-up tables)."""

import os

import numpy as np
import pytest
import xarray as xr
import zarr

EXPECTED_VARIABLES = ('Rrs', 'BRDFg', 'aot550', 'flags', 'mask', 'wind', 'sza', 'vza', 'raa')


@pytest.fixture(scope='module')
def l2a(l2a_zarr):
    with xr.open_zarr(l2a_zarr, group='0') as ds:
        yield ds.load()


@pytest.fixture(scope='module')
def l1c(l1c_path):
    with xr.open_dataset(l1c_path) as ds:
        yield ds.load()


def test_groups(l2a_zarr):
    groups = {name for name, _ in zarr.open_group(l2a_zarr, mode='r').groups()}
    assert {'0', 'ancillary'} <= groups


def test_zarr_encoding(l2a_zarr):
    """Every image variable is compressed (Blosc), the wind is packed like the other variables."""
    level = zarr.open_group(l2a_zarr, mode='r')['0']
    for var in EXPECTED_VARIABLES:
        codecs = [c.to_dict()['name'] for c in level[var].metadata.codecs]
        if codecs[0] == 'sharding_indexed':
            codecs = [c['name'] for c in level[var].metadata.codecs[0].to_dict()['configuration']['codecs']]
        assert 'blosc' in codecs, (var, codecs)
    assert level['wind'].dtype == np.int16
    assert level['wind'].attrs['units'] == 'm s-1'


def test_quicklooks(l2a_zarr, l2a):
    """PNG quicklooks next to the store: image size, transparent outside the image, water in colour."""
    from PIL import Image
    base = os.path.splitext(l2a_zarr)[0]
    images = {}
    for kind in ('l1c_rgb', 'l2a_rgb', 'sunglint'):
        with Image.open(f'{base}_ql_{kind}.png') as im:
            assert im.mode == 'RGBA'
            images[kind] = np.asarray(im)
        assert images[kind].shape == (l2a.sizes['y'], l2a.sizes['x'], 4)    # 128 px: no reduction
        assert (images[kind][..., 3] == 255).mean() > 0.9, kind
    water = (l2a.mask == 0).values
    rgb = images['l2a_rgb'][..., :3].astype(int)
    assert (rgb[water].std(axis=1) > 0).mean() > 0.5          # coloured water, not the grey of the land
    assert (rgb[~water] == round(0.6 * 255)).all()
    assert images['l1c_rgb'][..., :3].max() > 100


def test_variables_and_grid(l2a, l1c):
    missing = [v for v in EXPECTED_VARIABLES if v not in l2a]
    assert not missing
    np.testing.assert_allclose(l2a.x, l1c.x)
    np.testing.assert_allclose(l2a.y, l1c.y)
    assert l2a.Rrs.dims == ('wl', 'y', 'x')
    assert l2a.sizes['wl'] == 11


def test_attributes(l2a):
    assert l2a.attrs['aot_resolution'] == 240
    assert l2a.attrs['opac_model'] in ('ARCT_rh70', 'COAV_rh70', 'DESE_rh70', 'MACL_rh70', 'URBA_rh70', 'WASO_rh0')


def test_mask_and_flags(l2a):
    assert set(np.unique(l2a.mask.values)) <= {0, 1}
    assert l2a.flags.dtype == np.uint32
    assert 0.2 < float((l2a.mask == 0).mean()) < 0.95      # water, land and shoreline in the subset


def test_rrs_on_water(l2a):
    water = (l2a.mask == 0).values
    for wl in (443, 560, 665, 865):
        rrs = l2a.Rrs.sel(wl=wl).values
        assert np.isfinite(rrs[water]).mean() > 0.99, wl
    rrs560 = l2a.Rrs.sel(wl=560).values[water]
    assert 0 < np.median(rrs560) < 0.1


def test_aot(l2a):
    aot = l2a.aot550.values
    assert np.all(np.isfinite(aot))
    assert np.all(aot >= 0)                                  # clipped at 0 after the smoothing
    assert float(aot.max()) < 1.5
