"""Unit tests for Kernel.aerosol_visible_chunk."""

from unittest.mock import Mock
import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def kernel_visible():
    """Builds a minimal kernel to test aerosol_visible_chunk."""

    kernel = Kernel.__new__(Kernel)

    # --------------------------------------------------
    # General parameters
    # --------------------------------------------------
    kernel.pressure_ref = 101320.0
    kernel.aot_ref_lim = 0.8
    kernel.aot_ref_max = 0.2
    kernel.neg_pix_max = 0
    kernel.monoview = True
    kernel._cache = {}

    # --------------------------------------------------
    # Geometry and grid
    # --------------------------------------------------
    x = np.array([100.0, 120.0])
    y = np.array([200.0, 220.0])

    # Seven visible/NIR bands and two SWIR bands.
    wl = np.array([
        443.0, 490.0, 560.0, 665.0,
        705.0, 740.0, 783.0, 865.0, 1610.0
    ])

    # Constant reflectance to simplify tests
    bands = np.full((len(wl), 2, 2), 0.1, dtype=np.float32)

    angles = np.full((2, 2), 30.0)

    raster = xr.Dataset(
        data_vars={
            "bands": (("wl", "y", "x"), bands),
            "sza": (("y", "x"), angles),
            "vza": (("y", "x"), angles),
            "raa": (("y", "x"), angles),
            "water_pixel_number": (
                ("y", "x"),
                np.ones((2, 2), dtype=np.int32),
            ),
        },
        coords={
            "wl": wl,
            "x": x,
            "y": y,
        },
    )

    # --------------------------------------------------
    # Parameters and atmospheric LUT
    # --------------------------------------------------
    kernel.szas = np.array([30.0])
    kernel.vzas = np.array([30.0])
    kernel.azis = np.array([150.0])

    aot_refs = np.array([0.0, 0.1, 0.2, 0.4, 0.8])

    # LUT of diffuse reflectance
    lut_values = np.zeros(
        (len(wl), len(aot_refs), 1, 1, 1),
        dtype=np.float32,
    )

    for i, aot in enumerate(aot_refs):
        lut_values[:, i, 0, 0, 0] = aot * 0.01

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

    # LUT of aot.
    aot_lut_values = np.zeros(
        (len(wl), len(aot_refs)),
        dtype=np.float32,
    )

    for i, aot in enumerate(aot_refs):
        aot_lut_values[:, i] = aot * 0.1

    kernel.aot_lut = xr.DataArray(
        aot_lut_values,
        dims=("wl", "aot_ref"),
        coords={
            "wl": wl,
            "aot_ref": aot_refs,
        },
    )

    # --------------------------------------------------
    # Gaseous transmittances
    # --------------------------------------------------
    kernel.Tg_diff_raster = xr.DataArray(
        np.ones((len(wl), 2, 2), dtype=np.float32),
        dims=("wl", "y", "x"),
        coords={"wl": wl, "x": x, "y": y},
    )

    kernel.Tg_raster = xr.DataArray(
        np.ones((len(wl), 2, 2), dtype=np.float32),
        dims=("wl", "y", "x"),
        coords={"wl": wl, "x": x, "y": y},
    )

    # --------------------------------------------------
    # Rayleigh and sunglint
    # --------------------------------------------------
    kernel.rot = xr.DataArray(
        np.full(len(wl), 0.01, dtype=np.float32),
        dims=("wl",),
        coords={"wl": wl},
    )

    kernel.sunglint_eps = xr.DataArray(
        np.ones(len(wl), dtype=np.float32),
        dims=("wl",),
        coords={"wl": wl},
    )

    kernel.wind_img = xr.Dataset(
        {
            "wind": (
                ("y", "x"),
                np.full((2, 2), 2.0, dtype=np.float32),
            )
        },
        coords={"x": x, "y": y},
    )

    # --------------------------------------------------
    # CAMS data
    # --------------------------------------------------
    kernel.cams = Mock()

    kernel.cams.raster.sp = xr.DataArray(
        np.full((2, 2), 101320.0, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
    )

    kernel.cams.cams_aod = xr.DataArray(
        np.full((3,2,2), 0.2, dtype=np.float32),
        dims=("wl", "y", "x"),
        coords={
            "wl": [440.0, 550.0, 670.0],
            "x": x,
            "y": y,
        },
    )

    # --------------------------------------------------
    # Simulated atmospheric interpolation
    # --------------------------------------------------
    kernel._R_ = Mock()

    # the result will be replaced in individual tests
    # to control the calculated reflectances
    kernel._R_.interp_Rlut_rayleigh.return_value = np.zeros(
        (len(wl), 2, 2),
        dtype=np.float32,
    )

    # Spectral multiplication used to build the sunglint
    def fake_multiplicate(spectrum, field, output):
        output[:] = spectrum[:, None, None] * field[None, :, :]
        return output

    kernel._R_._multiplicate.side_effect = fake_multiplicate

    kernel.test_raster = raster

    return kernel

# ============================================================
# Tests AOT Selection
# ============================================================

class TestAerosolVisibleChunkAOTSelection:
    
    def test_fallback_cams(self,kernel_visible, monkeypatch):
        """Use CAMS when no AOT candidate has been accepted"""
        kernel = kernel_visible
        raster = kernel_visible.test_raster
        
        monkeypatch.setattr("grs.grs_kernel._aot_ref_candidates", lambda *args, **kwargs: np.array([]))
        kernel.cams.cams_aod.loc[{"wl": 550.0}] = 0.35
        result = kernel_visible.aerosol_visible_chunk(raster)
        
        assert result is not None
        assert isinstance(result, xr.Dataset)

        #AOD Cams at 550nm is equal to 0.2 on each pixel
        #its spatial mean must be equal to 0.2 as well
        assert np.isclose(
            result["aot_ref"].values[0,0],
            0.35
        )
        
        assert result["number_of_pixel"].values[0,0] == 4
        
    def test_selects_larges_valid_aot(self,kernel_visible, monkeypatch):
        """Selects the largest AOT candidate amongst admissible candidates"""
        kernel = kernel_visible
        raster = kernel_visible.test_raster
        
        expected_candidates = np.array([0.0,0.1,0.3])
        monkeypatch.setattr("grs.grs_kernel._aot_ref_candidates", lambda *args, **kwargs: expected_candidates)
        
        result = kernel_visible.aerosol_visible_chunk(raster)
        
        assert np.isclose(
            result["aot_ref"].values[0,0],
            0.3
        )
        
        assert result["number_of_pixel"].values[0,0] == 4    

# ============================================================
# Tests Outputs and edge cases
# ============================================================

class TestAerosolVisibleChunkOutputs:
        
    def test_empty_block_returns_None(self,kernel_visible):
        raster = kernel_visible.test_raster.isel(x=slice(0, 0))
        
        result = kernel_visible.aerosol_visible_chunk(raster)

        assert result is None

    def test_nominal(self,kernel_visible):
        """Checks that aerosol_visible_chunk returns a valid output"""

        kernel = kernel_visible
        raster = kernel.test_raster

        result = kernel.aerosol_visible_chunk(raster)

        assert result is not None
        assert isinstance(result, xr.Dataset)

        assert "aot_ref" in result
        assert "number_of_pixel" in result

        assert result["aot_ref"].shape == (1, 1)
        assert result["number_of_pixel"].shape == (1, 1)

        assert np.isfinite(result["aot_ref"].values).all()
        assert result["number_of_pixel"].values[0, 0] == 4

    def test_keeps_xres(self,kernel_visible, monkeypatch):
        """Stores reflectance arrays for every tested AOT when keep_xres=True."""
        kernel = kernel_visible

        monkeypatch.setattr(
            "grs.grs_kernel._aot_ref_candidates",
            lambda *args, **kwargs: np.array([0.1]),
        )

        kernel.aerosol_visible_chunk(
            kernel.test_raster,
            keep_xres=True,
        )

        assert hasattr(kernel, "xres")
        assert "Rrs" in kernel.xres
        assert kernel.xres.sizes["aot_ref"] == 21
        assert kernel.xres.sizes["wl"] == 9
        assert kernel.xres.sizes["y"] == 2
        assert kernel.xres.sizes["x"] == 2