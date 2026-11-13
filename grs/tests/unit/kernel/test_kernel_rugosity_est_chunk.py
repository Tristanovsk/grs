"""Unit tests for Kernel.rugosity_est_chunk."""

import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel

# ============================================================
# Fixtures and helpers
# ============================================================

class FakeRasterization:
    """Replace the interpolation of the atmospheric LUT"""

    def __init__(self, diffuse_value=0.0):
        self.diffuse_value = diffuse_value
        self.calls = []

    def interp_Rlut_rayleigh(
        self,
        szas,
        sza,
        vzas,
        vza,
        azis,
        azi,
        nwl,
        ny,
        nx,
        lut_values,
    ):
        self.calls.append({
            "sza": sza.copy(),
            "vza": vza.copy(),
            "azi": azi.copy(),
            "nwl": nwl,
            "ny": ny,
            "nx": nx,
        })

        return np.full(
            (nwl, ny, nx),
            self.diffuse_value,
            dtype=np.float64,
        )

def make_test_raster():
        """Create an artificial raster with 4 bands"""

        wl = np.array([700.0, 800.0, 900.0, 1000.0])
        y = np.array([10.0, 20.0])
        x = np.array([100.0, 110.0, 120.0])

        bands = np.array([
            np.full((2, 3), 0.10),
            np.full((2, 3), 0.12),
            np.full((2, 3), 0.14),
            np.full((2, 3), 0.16),
        ])

        sza = np.full((4, 2, 3), 30.0)
        vza = np.full((4, 2, 3), 10.0)
        raa = np.full((4, 2, 3), 180.0)

        return xr.Dataset(
            data_vars={
                "bands": (("wl", "y", "x"), bands),
                "sza": (("wl", "y", "x"), sza),
                "vza": (("wl", "y", "x"), vza),
                "raa": (("wl", "y", "x"), raa),
            },
            coords={
                "wl": wl,
                "y": y,
                "x": x,
            },
        )
        
@pytest.fixture
def kernel_and_raster():
    """Create a minimal test Kernel and raster"""

    kernel = Kernel.__new__(Kernel)

    raster = make_test_raster()

    kernel.monoview = False
    kernel.pressure_ref = 101320.0

    kernel.wl_true = raster.wl

    kernel.szas = np.array([30.0])
    kernel.vzas = np.array([10.0])
    kernel.azis = np.array([0.0])

    # Artificial LUT diffuse.
    kernel.Rdiff_lut = xr.DataArray(
        np.zeros((4, 2)),
        dims=("wl", "aot_ref"),
        coords={
            "wl": raster.wl,
            "aot_ref": [0.0, 0.8],
        },
    )

    pressure = xr.DataArray(
        np.full((2, 3), 101320.0),
        dims=("y", "x"),
        coords={
            "y": raster.y,
            "x": raster.x,
        },
    )

    kernel._cache = {"sp": pressure}

    transmittance = xr.DataArray(
        np.ones((4, 2, 3)),
        dims=("wl", "y", "x"),
        coords={
            "wl": raster.wl,
            "y": raster.y,
            "x": raster.x,
        },
    )

    kernel.Tg_raster = transmittance
    kernel.Tg_diff_raster = transmittance

    kernel._R_ = FakeRasterization(diffuse_value=0.0)

    return kernel, raster

# ============================================================
# Tests
# ============================================================

class TestRugosityEstChunk:

    def test_output_structure(self,kernel_and_raster):
        kernel, raster = kernel_and_raster

        result = kernel.rugosity_est_chunk(
            raster,
            winds=np.array([1.0, 5.0, 10.0]),
            idx_nir=[0, 1, 2, 3],
        )

        assert isinstance(result, xr.Dataset)

        assert set(result.data_vars) == {
            "wind",
            "cost",
            "brdf_g",
        }

        for variable in ("wind", "cost", "brdf_g"):
            assert result[variable].dims == ("y", "x")
            assert result[variable].shape == (2, 3)

        np.testing.assert_array_equal(
            result.x.values,
            raster.x.values,
        )

        np.testing.assert_array_equal(
            result.y.values,
            raster.y.values,
        )

        assert result.attrs["units"] == "m.s-1"

    def test_selects_first_wind_when_costs_are_equal_when_p_slope_is_wind_independant(self,kernel_and_raster):
        # This test will need to be readapted if p_slope wind dependance is reactivated !!!
        
        kernel, raster = kernel_and_raster

        winds = np.array([2.0, 5.0, 8.0])

        result = kernel.rugosity_est_chunk(
            raster,
            winds=winds,
            idx_nir=[0, 1, 2, 3],
        )

        np.testing.assert_allclose(
            result["wind"].values,
            np.full((2, 3), winds[0]),
        )

        np.testing.assert_allclose(
            result["cost"].values,
            0.0,
            atol=1e-12,
        )

    def test_subtracts_diffuse_reflectance(self,kernel_and_raster):
        kernel, raster = kernel_and_raster

        kernel._R_ = FakeRasterization(diffuse_value=0.01)

        result = kernel.rugosity_est_chunk(
            raster,
            winds=np.array([2.0, 5.0]),
            idx_nir=[0, 1, 2, 3],
        )

        expected_brdf_g = 0.16 - 0.01

        np.testing.assert_allclose(
            result["brdf_g"].values,
            expected_brdf_g,
        )

    def test_azimuth_conversion(self,kernel_and_raster):
        kernel, raster = kernel_and_raster
        
        kernel.rugosity_est_chunk(
            raster,
            winds=np.array([2.0]),
            idx_nir=[0,1,2,3],
        )
        
        assert len(kernel._R_.calls) == 1
        call = kernel._R_.calls[0]
        
        np.testing.assert_allclose(
            call["azi"],
            0,0
        )
        
        assert call["nwl"] == 4
        assert call["ny"] == 2
        assert call["nx"] == 3

    @pytest.mark.parametrize(
        "empty_dimension",
        ["x", "y"],
    )
    def test_empty_block(self,kernel_and_raster, empty_dimension):
        kernel, raster = kernel_and_raster
        if empty_dimension == "x":
            empty_raster = raster.isel(x=slice(0,0))
        else:
            empty_raster = raster.isel(y=slice(0,0))
        
        result = kernel.rugosity_est_chunk(
            empty_raster,
            winds=np.array([2.0,5.0]),
            idx_nir=[0,1,2,3]
        )
        
        assert result is None

    def test_nan_angles(self,kernel_and_raster):
        kernel, raster = kernel_and_raster
        
        #Invalid angles at pixel (0,0)
        raster["sza"].values[:,0,0] = np.nan
        result = kernel.rugosity_est_chunk(
            raster,
            winds=np.array([2.0,5.0]),
            idx_nir=[0,1,2,3]
        )
        
        assert np.isnan(result["wind"].values[0,0])
        assert np.isnan(result["cost"].values[0,0])
        
        #Other pixels have a valid result
        assert np.isfinite(result["wind"].values[1,1])
        assert np.isfinite(result["cost"].values[1,1])