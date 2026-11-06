"""Unit tests for Kernel.get_coarse_masked_raster."""

import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def kernel_coarse():
    """Create a minimal Kernel instance for testing coarse raster generation."""
    instance = Kernel.__new__(Kernel)

    instance.raster = xr.Dataset(
        data_vars={
            "sza": (
                ("y", "x"),
                [
                    [10., 20., 30., 40.],
                    [30., 40., 50., 60.],
                    [50., 60., 70., 80.],
                    [70., 80., 90., 100.],
                ],
            ),
            "vza": (
                ("y", "x"),
                [
                    [1., 2., 3., 4.],
                    [3., 4., 5., 6.],
                    [5., 6., 7., 8.],
                    [7., 8., 9., 10.],
                ],
            ),
            "raa": (
                ("y", "x"),
                [
                    [100., 110., 120., 130.],
                    [110., 120., 130., 140.],
                    [120., 130., 140., 150.],
                    [130., 140., 150., 160.],
                ],
            ),
            "bands": (
                ("y", "x"),
                [
                    [0.1, 0.2, 0.3, 0.4],
                    [0.3, 0.4, 0.5, 0.6],
                    [0.5, 0.6, 0.7, 0.8],
                    [0.7, 0.8, 0.9, 1.0],
                ],
            ),
        },
        coords={
            "x": [10, 30, 50, 70],
            "y": [10, 30, 50, 70],
        },
    )

    return instance

@pytest.fixture
def kernel_coarse_odd():
    """Create a 3x3 raster to test incomplete coarse blocks."""
    instance = Kernel.__new__(Kernel)

    instance.raster = xr.Dataset(
        data_vars={
            "sza": (
                ("y", "x"),
                [
                    [10., 20., 30.],
                    [30., 40., 50.],
                    [50., 60., 70.],
                ],
            ),
            "vza": (
                ("y", "x"),
                np.ones((3, 3)),
            ),
            "raa": (
                ("y", "x"),
                np.ones((3, 3)),
            ),
            "bands": (
                ("y", "x"),
                np.ones((3, 3)),
            ),
        },
        coords={
            "x": [10, 30, 50],
            "y": [10, 30, 50],
        },
    )

    return instance

# ============================================================
# Tests Aggregation
# ============================================================

class TestCoarseMaskedRasterAggregation:
    
    def test_mean(self,kernel_coarse):
        """Test that raster variables are averaged over coarse pixel blocks."""

        result = kernel_coarse.get_coarse_masked_raster(
            resolution=40,
            min_pixel_prop=0
        )

        expected_sza = np.array([
            [25., 45.],
            [65., 85.],
        ])

        np.testing.assert_allclose(
            result["sza"].values,
            expected_sza
        )

        assert result.sizes["x"] == 2
        assert result.sizes["y"] == 2

    def test_incomplete_blocks(self,kernel_coarse_odd):
        """Test that incomplete blocks at the raster edges are padded."""

        result = kernel_coarse_odd.get_coarse_masked_raster(
            resolution=40,
            min_pixel_prop=0
        )

        expected_sza = np.array([
            [25., 40.],
            [55., 70.],
        ])

        np.testing.assert_allclose(
            result["sza"].values,
            expected_sza
        )

        expected_counts = np.array([
            [4, 2],
            [2, 1],
        ])

        np.testing.assert_array_equal(
            result["water_pixel_number"].values,
            expected_counts
        )
        
# ============================================================
# Tests Masking and pixel counting
# ============================================================

class TestCoarseMaskedRasterMasking:
    
    def test_with_mask(self,kernel_coarse):
        """Test that masked pixels are excluded from coarse averages."""

        mask = xr.DataArray(
            [
                [0, 1, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ],
            dims=("y", "x"),
            coords=kernel_coarse.raster.coords,
        )

        result = kernel_coarse.get_coarse_masked_raster(
            mask=mask,
            resolution=40,
            min_pixel_prop=0
        )

        # The pixel with sza=20 is masked.
        # The remaining values are 10, 30 and 40.
        expected_mean = (10. + 30. + 40.) / 3.

        np.testing.assert_allclose(
            result["sza"].values[0, 0],
            expected_mean
        )

        assert result["water_pixel_number"].values[0, 0] == 3

    def test_water_pixel_number(self,kernel_coarse):
        """Test the number of valid pixels in each coarse block."""

        result = kernel_coarse.get_coarse_masked_raster(
            resolution=40,
            min_pixel_prop=0
        )

        expected = np.full((2, 2), 4)

        np.testing.assert_array_equal(
            result["water_pixel_number"].values,
            expected
        )
        
    def test_min_pixel_prop(self,kernel_coarse):
        """Test filtering of coarse blocks according to their valid pixel count."""

        mask = xr.DataArray(
            [
                [0, 0, 0, 1],
                [0, 0, 1, 1],
                [0, 1, 0, 1],
                [1, 1, 1, 0],
            ],
            dims=("y", "x"),
            coords=kernel_coarse.raster.coords,
        )

        result = kernel_coarse.get_coarse_masked_raster(
            mask=mask,
            resolution=40,
            min_pixel_prop=0.5
        )

        expected_counts = np.array([
            [4, 1],
            [1, 2],
        ])

        np.testing.assert_array_equal(
            kernel_coarse.coarse_masked_raster["water_pixel_number"].values,
            expected_counts
        )

        assert result["sza"].values[0, 0] == pytest.approx(25.)

        assert np.isnan(result["sza"].values[0, 1])
        assert np.isnan(result["sza"].values[1, 0])
        assert np.isnan(result["sza"].values[1, 1])

# ============================================================
# Tests Output and attributes
# ============================================================

class TestCoarseMaskedRasterOutputs:
    
    def test_attributes(self,kernel_coarse):
        """Test that the coarse raster attributes are correctly updated."""

        result = kernel_coarse.get_coarse_masked_raster(
            resolution=40,
            min_pixel_prop=0
        )

        assert kernel_coarse.coarse_resolution == 40
        assert kernel_coarse.coarse_masked_raster.attrs["resolution"] == 40
        assert result.attrs["resolution"] == 40


