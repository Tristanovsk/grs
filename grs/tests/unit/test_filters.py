"""Smoothing filters of the AOT retrieval (Kernel.filter2d, Kernel.smoothing, Kernel.conv_mapping)."""

import numpy as np
import pytest
import xarray as xr
from scipy import ndimage

from grs.grs_kernel import Kernel


def test_filter2d_weighted_mean():
    image = np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])
    weight = np.ones_like(image)
    weight[0, 0] = 0                       # null weight: ignored by the neighbours, value kept
    out = Kernel.filter2d(image, weight, np.array([3, 3]))
    assert out[0, 0] == 1.
    assert out[1, 1] == pytest.approx(np.mean(image.ravel()[1:]))
    # window truncated at the edges
    assert out[2, 2] == pytest.approx(np.mean([5., 6., 8., 9.]))
    np.testing.assert_array_equal(image, np.arange(1., 10.).reshape(3, 3))   # input unchanged


def test_filter2d_weights_and_nan():
    image = np.array([[1., 100.], [3., 5.]])
    weight = np.array([[1., np.nan], [3., 1.]])
    out = Kernel.filter2d(image, weight, np.array([3, 3]))
    assert out[0, 1] == 100.               # NaN weight: kept, and ignored by the neighbours
    assert out[0, 0] == pytest.approx((1 * 1 + 3 * 3 + 1 * 5) / 5)


def test_filter2d_one_block_window_is_identity():
    image = np.random.default_rng(0).random((4, 5))
    np.testing.assert_array_equal(Kernel.filter2d(image, np.ones_like(image), np.array([1, 1])), image)


def test_smoothing_box_equals_conv_mapping():
    rng = np.random.default_rng(1)
    values = rng.random((20, 25))
    values[rng.random(values.shape) < 0.2] = np.nan
    raster = xr.Dataset({'aot_ref': (['y', 'x'], values)}, coords=dict(y=np.arange(20.), x=np.arange(25.)))
    kernel = Kernel.__new__(Kernel)        # the smoothing only uses static methods
    box = kernel.smoothing(raster, varname='aot_ref', mask=np.ones((5, 5)), windows=np.array([5, 5]))
    weights = 1 + 0 * (1 / values ** 2)
    first = Kernel.filter2d(values, weights, np.array([5, 5]))
    ref = ndimage.generic_filter(first, function=Kernel.conv_mapping, footprint=np.ones((5, 5)), mode='nearest')
    np.testing.assert_allclose(box.aot_ref.values, ref, rtol=0, atol=1e-12)
