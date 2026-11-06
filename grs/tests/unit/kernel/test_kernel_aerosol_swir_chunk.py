"""Unit tests for Kernel.aerosol_swir_chunk."""

from unittest.mock import Mock
import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel


# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def kernel_swir():
    """Builds a minimal Kernel to test aerosol_swir_chunk"""

    kernel = Kernel.__new__(Kernel)

    kernel.pressure_ref = 101320.0
    kernel.aot_ref_lim = 0.8
    kernel.aot_ref_max = 0.5
    kernel._cache = {}

    x = np.array([100.0, 120.0])
    y = np.array([200.0, 220.0])
    wl = np.array([560.0, 865.0, 1610.0])

    bands = np.array([
        [[0.50, 0.50], [0.50, 0.50]],
        [[0.10, 0.10], [0.10, 0.10]],
        [[0.08, 0.08], [0.08, 0.08]],
    ])

    angles = np.full((2, 2), 30.0)

    raster = xr.Dataset(
        data_vars={
            "bands": (("wl", "y", "x"), bands),
            "sza": (("y", "x"), angles),
            "vza": (("y", "x"), angles),
            "raa": (("y", "x"), angles),
        },
        coords={"wl": wl, "x": x, "y": y},
    )

    kernel.cams = Mock()
    kernel.cams.raster.sp = xr.DataArray(
        np.full((2, 2), 101320.0),
        dims=("y", "x"),
        coords={"x": x, "y": y},
    )

    kernel._R_ = Mock()

    kernel.szas = np.array([30.0])
    kernel.vzas = np.array([30.0])
    kernel.azis = np.array([150.0])
    
    aot_refs = np.array([0.0, 0.1, 0.2])
    lut_values = np.zeros((3, 3, 1, 1, 1))

    for i, aot in enumerate(aot_refs):
        lut_values[:, i, 0, 0, 0] = aot

    kernel.Rdiff_lut = xr.DataArray(
        lut_values,
        dims=("wl", "aot_ref", "sza", "vza", "azi"),
        coords={
            "wl": wl,
            "aot_ref": aot_refs,
            "sza": [30.0],
            "vza": [30.0],
            "azi": [150.0],
        },
    )

    kernel.Tg_diff_raster = xr.DataArray(
        np.ones((3, 2, 2)),
        dims=("wl", "y", "x"),
        coords={"wl": wl, "x": x, "y": y},
    )

    kernel.test_raster = raster

    return kernel

def configure_diffuse_interpolation(kernel, diffuse_by_aot):
    """Configures the mock to send a known diffuse contribution"""

    def fake_interp(*args):
        lut_values = np.asarray(args[-1])
        aot = float(lut_values.flat[0])

        diffuse = diffuse_by_aot[aot]

        return np.full((2, 2, 2), diffuse, dtype=float)

    kernel._R_.interp_Rlut_rayleigh.side_effect = fake_interp


# ============================================================
# Tests AOT selection and fallback
# ============================================================

class TestAerosolSwirChunkAOTSelection:
    
    def test_selects_largest_valid_aot(self,kernel_swir):
        configure_diffuse_interpolation(
            kernel_swir,
            {
                0.0: 0.0,
                0.1: 0.02,
                0.2: 0.09,
            },
        )

        aot_refs = [0.0, 0.1, 0.2]

        result = kernel_swir.aerosol_swir_chunk(
            kernel_swir.test_raster,
            aot_refs=aot_refs,
            idx_swir=[-2, -1],
            neg_pix_rel_thresh=0.05,
        )

        assert result == pytest.approx(0.1)
        assert kernel_swir.aot_ref_max == pytest.approx(0.1)

    def test_uses_fallback_when_all_counts_are_zero(self,kernel_swir):
        configure_diffuse_interpolation(
            kernel_swir,
            {
                0.0: 0.0,
                0.1: 0.0,
                0.2: 0.0,
            },
        )

        result = kernel_swir.aerosol_swir_chunk(
            kernel_swir.test_raster,
            aot_refs=[0.0, 0.1, 0.2],
        )

        assert result == pytest.approx(0.16)
        assert kernel_swir.aot_ref_max == pytest.approx(0.16)


    def test_selects_largest_valid_aot(self,kernel_swir):
        configure_diffuse_interpolation(
            kernel_swir,
            {
                0.0: 0.0,
                0.1: 0.02,
                0.2: 0.09,
            },
        )

        aot_refs = [0.0, 0.1, 0.2]

        result = kernel_swir.aerosol_swir_chunk(
            kernel_swir.test_raster,
            aot_refs=aot_refs,
            idx_swir=[-2, -1],
            neg_pix_rel_thresh=0.05,
        )

        assert result == pytest.approx(0.1)
        assert kernel_swir.aot_ref_max == pytest.approx(0.1)
    
    def test_enforces_lut_limit(self,kernel_swir):
            kernel_swir.aot_ref_lim = 0.2
            kernel_swir.aot_ref_cams_max = 0.15
    
            configure_diffuse_interpolation(
                kernel_swir,
                {
                    0.0: 0.0,
                    0.1: 0.02,
                    0.2: 0.09,
                },
            )
    
            result = kernel_swir.aerosol_swir_chunk(
                kernel_swir.test_raster,
                aot_refs=[0.0, 0.1, 0.2],
                neg_pix_rel_thresh=2.0,
            )
    
            assert result == pytest.approx(0.15)
            assert kernel_swir.aot_ref_max == pytest.approx(0.15)

# ============================================================
# Tests Pixel filtering
# ============================================================

class TestAerosolSwirChunkPixelFiltering:

    def test_ignores_pixels_negative_at_zero_aot(self,kernel_swir):
        diffuse_by_aot = {
            0.0: np.array([
                [[0.00, 0.00], [0.00, 0.00]],
                [[0.10, 0.00], [0.00, 0.00]],
            ]),
            0.1: np.array([
                [[0.00, 0.02], [0.00, 0.00]],
                [[0.10, 0.02], [0.00, 0.00]],
            ]),
            0.2: np.array([
                [[0.09, 0.09], [0.09, 0.09]],
                [[0.09, 0.09], [0.09, 0.09]],
            ]),
        }

        def fake_interp(*args):
            aot = float(np.asarray(args[-1]).flat[0])
            return diffuse_by_aot[aot]

        kernel_swir._R_.interp_Rlut_rayleigh.side_effect = fake_interp

        result = kernel_swir.aerosol_swir_chunk(
            kernel_swir.test_raster,
            aot_refs=[0.0, 0.1, 0.2],
            neg_pix_rel_thresh=0.5,
        )

        assert result == pytest.approx(0.1)


# ============================================================
# Tests Ouptuts
# ============================================================

class TestAerosolSwirChunkOutputs:

    def test_empty_raster(self,kernel_swir):
        raster = kernel_swir.test_raster.isel(x=slice(0, 0))

        result = kernel_swir.aerosol_swir_chunk(raster)

        assert result is None
        kernel_swir._R_.interp_Rlut_rayleigh.assert_not_called()

    def test_keep_xres(self,kernel_swir):
        configure_diffuse_interpolation(
            kernel_swir,
            {
                0.0: 0.0,
                0.1: 0.02,
                0.2: 0.09,
            },
        )

        aot_refs = [0.0, 0.1, 0.2]

        kernel_swir.aerosol_swir_chunk(
            kernel_swir.test_raster,
            aot_refs=aot_refs,
            keep_xres=True,
        )

        assert hasattr(kernel_swir, "xres")

        assert "Rrs" in kernel_swir.xres
        assert "aot_ref" in kernel_swir.xres.dims

        assert kernel_swir.xres.sizes["aot_ref"] == 3
        assert kernel_swir.xres.sizes["wl"] == 2
        assert kernel_swir.xres.sizes["y"] == 2
        assert kernel_swir.xres.sizes["x"] == 2

        np.testing.assert_allclose(
            kernel_swir.xres.aot_ref.values,
            aot_refs,
        )

    def test_keep_xres(self,kernel_swir):
        configure_diffuse_interpolation(
            kernel_swir,
            {
                0.0: 0.0,
                0.1: 0.02,
                0.2: 0.09,
            },
        )

        aot_refs = [0.0, 0.1, 0.2]

        kernel_swir.aerosol_swir_chunk(
            kernel_swir.test_raster,
            aot_refs=aot_refs,
            keep_xres=True,
        )

        assert hasattr(kernel_swir, "xres")

        assert "Rrs" in kernel_swir.xres
        assert "aot_ref" in kernel_swir.xres.dims

        assert kernel_swir.xres.sizes["aot_ref"] == 3
        assert kernel_swir.xres.sizes["wl"] == 2
        assert kernel_swir.xres.sizes["y"] == 2
        assert kernel_swir.xres.sizes["x"] == 2

        np.testing.assert_allclose(
            kernel_swir.xres.aot_ref.values,
            aot_refs,
        )

