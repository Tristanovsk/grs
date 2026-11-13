"""Unit tests for Kernel.coarsening_factor."""

import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def kernel():
    """Create a Kernel instance without calling its constructor."""
    instance = Kernel.__new__(Kernel)

    instance.raster = xr.Dataset(
        coords={
            "x": [10, 30, 50, 70, 90],
            "y": [100, 120, 140, 160, 180],
        }
    )

    return instance

# ============================================================
# Test coarsening factor
# ============================================================

class TestCoarseningFactor:
    
    def test_standard_resolution(self,kernel):
        """Test coarsening with a target resolution that is a multiple
        of the pixel size.
        """
        result = kernel.coarsening_factor(resolution=60)

        assert result == {"x": 3, "y": 3}


    def test_different_pixel_sizes(self):
        """Test coarsening when x and y have different pixel sizes."""
        instance = Kernel.__new__(Kernel)

        instance.raster = xr.Dataset(
            coords={
                "x": [0, 20, 40, 60, 80],
                "y": [0, 30, 60, 90, 120],
            }
        )

        result = instance.coarsening_factor(resolution=120)

        assert result == {"x": 6, "y": 4}
        
    @pytest.mark.parametrize(
        "resolution, expected",
        [
            (20, {"x": 1, "y": 1}),
            (40, {"x": 2, "y": 2}),
            (60, {"x": 3, "y": 3}),
            (120, {"x": 6, "y": 6}),
        ],
    )
    def test_multiple_resolutions(self,kernel, resolution, expected):
        """Test several valid target resolutions."""
        result = kernel.coarsening_factor(resolution)

        assert result == expected


    @pytest.mark.parametrize("resolution", [10, 30, 50, 70])
    def test_invalid_resolution(self,kernel, resolution):
        """Test that a resolution that is not a multiple of the pixel size
        raises a ValueError.
        """
        with pytest.raises(ValueError, match="not a multiple"):
            kernel.coarsening_factor(resolution)


    @pytest.mark.parametrize("resolution", [0, -20])
    def test_non_positive_resolution(self,kernel, resolution):
        """Test that a non-positive resolution raises a ValueError."""
        with pytest.raises(ValueError, match="not a multiple"):
            kernel.coarsening_factor(resolution)


    def test_descending_coordinates(self):
        """Test that descending coordinates are handled correctly."""
        instance = Kernel.__new__(Kernel)

        instance.raster = xr.Dataset(
            coords={
                "x": [90, 70, 50, 30, 10],
                "y": [180, 160, 140, 120, 100],
            }
        )

        result = instance.coarsening_factor(resolution=60)

        assert result == {"x": 3, "y": 3}