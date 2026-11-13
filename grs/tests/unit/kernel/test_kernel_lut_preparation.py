"""Unit tests for kernel.lut_preparation"""


import numpy as np
import pytest
import xarray as xr
import logging

from types import SimpleNamespace

from grs.grs_kernel import Kernel
import grs.grs_kernel as kernel_module

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def mock_auxdata(monkeypatch):
    """Replace AuxData with synthetic auxiliary data."""

    def fake_auxdata(wl=None):
        if wl is None:
            wl = np.array([500., 600.])

        wl = np.asarray(wl)

        return SimpleNamespace(
            sunglint_eps=xr.DataArray(
                np.ones(len(wl)),
                dims=("wl",),
                coords={"wl": wl},
            ),
            rot=xr.DataArray(
                np.full(len(wl), 0.1),
                dims=("wl",),
                coords={"wl": wl},
            ),
        )

    monkeypatch.setattr(kernel_module, "AuxData", fake_auxdata)

@pytest.fixture
def kernel_lut():
    """Create a Kernel with synthetic aerosol LUTs."""

    wl_true = np.array([500., 600., 700., 800.])

    raster = xr.Dataset(
        data_vars={
            "sza": (
                ("y", "x"),
                [
                    [30., 30.],
                    [40., 40.],
                ],
            ),
            "vza": (
                ("y", "x"),
                [
                    [10., 10.],
                    [20., 20.],
                ],
            ),
            "raa": (
                ("y", "x"),
                [
                    [0., 0.],
                    [90., 90.],
                ],
            ),
            "bands": (
                ("y", "x", "wl"),
                np.ones((2, 2, 4)),
            ),
        },
        coords={
            "x": [0, 20],
            "y": [0, 20],
            "wl": wl_true,
        },
    )

    raster = raster.assign_coords(wl_true=("wl", wl_true))

    prod = SimpleNamespace(
        raster=raster,
        sensor="S2",
        _type=np.float32,
    )

    winds = [1., 2., 3.]
    models = kernel_module.AEROSOL_MODELS
    aot_refs = np.linspace(0., 0.8, 5)
    szas = [0., 10., 20., 30., 40., 50.]
    vzas = [0., 10., 20., 30.]
    azis = [0., 90., 180., 270.]
    wls = [400., 500., 600., 700., 800., 900.]

    dims = (
        "wind",
        "model",
        "aot_ref",
        "sza",
        "vza",
        "azi",
        "wl",
    )

    coords = {
        "wind": winds,
        "model": models,
        "aot_ref": aot_refs,
        "sza": szas,
        "vza": vzas,
        "azi": azis,
        "wl": wls,
    }

    shape = tuple(len(coords[dim]) for dim in dims)

    I = np.full(shape, 0.2, dtype=float)

    aot = np.ones(shape, dtype=float)

    aero_lut = xr.Dataset(
        data_vars={
            "I": (dims, I),
            "aot": (dims, aot),
        },
        coords=coords,
    )

    trans_lut = xr.Dataset(
        data_vars={
            "Ttot_Ed": (dims, np.ones(shape, dtype=float)),
        },
        coords=coords,
    )

    kernel = Kernel.__new__(Kernel)

    kernel.prod = prod
    kernel.raster = raster
    kernel.aero_lut = aero_lut
    kernel.trans_lut = trans_lut
    kernel.wl_true = wl_true
    kernel._cache = {}

    return kernel

@pytest.fixture
def kernel_lut_spectral(kernel_lut):
    """Create a LUT with a wavelength-dependent intensity."""

    kernel = kernel_lut

    lut = kernel.aero_lut["I"]
    wls = lut["wl"].values

    intensities = np.broadcast_to(
        (wls / 1000)[None, None, None, None, None, None, :],
        lut.shape,
    ).copy()

    kernel.aero_lut["I"] = (lut.dims, intensities)

    return kernel

@pytest.fixture
def kernel_lut_spectral(kernel_lut):
    """Create a LUT with a wavelength-dependent intensity."""

    kernel = kernel_lut

    lut = kernel.aero_lut["I"]
    wls = lut["wl"].values

    intensities = np.broadcast_to(
        (wls / 1000)[None, None, None, None, None, None, :],
        lut.shape,
    ).copy()

    kernel.aero_lut["I"] = (lut.dims, intensities)

    return kernel

@pytest.fixture
def kernel_lut_diffuse(kernel_lut):
    """Create a Kernel with distinct intensities for each aerosol model."""

    kernel = kernel_lut

    dims = kernel.aero_lut["I"].dims
    shape = kernel.aero_lut["I"].shape
    intensities = np.broadcast_to(
        np.arange(1, 7, dtype=float)[None, :, None, None, None, None, None],
        shape,
    ).copy()

    kernel.aero_lut["I"] = (dims, intensities)

    return kernel

@pytest.fixture
def kernel_lut_aot(kernel_lut):
    """Create a Kernel with synthetic AOT values varying by model,
    reference AOT and wavelength.
    """

    kernel = kernel_lut

    lut = kernel.aero_lut["aot"]

    model_factors = np.arange(1, 7, dtype=float).reshape(
        1, 6, 1, 1, 1, 1, 1
    )

    aot_refs = lut.aot_ref.values
    wls = lut.wl.values

    base = (
        aot_refs[None, None, :, None, None, None, None]
        + wls[None, None, None, None, None, None, :] / 1000
    )

    values = np.broadcast_to(
        model_factors * base,
        lut.shape,
    ).copy()

    kernel.aero_lut["aot"] = (lut.dims, values)

    return kernel

@pytest.fixture
def kernel_lut_transmittance(kernel_lut):
    """Create a transmittance LUT with distinct values per aerosol model."""

    kernel = kernel_lut

    dims = kernel.trans_lut["Ttot_Ed"].dims
    shape = kernel.trans_lut["Ttot_Ed"].shape

    transmittances = np.broadcast_to(
        np.arange(1, 7, dtype=float)[None, :, None, None, None, None, None],
        shape,
    ).copy()

    kernel.trans_lut["Ttot_Ed"] = (dims, transmittances)
    return kernel

@pytest.fixture
def kernel_lut_wind_dependent(kernel_lut):
    """Create an aerosol LUT whose intensity depends on wind speed."""

    kernel = kernel_lut

    dims = kernel.aero_lut["I"].dims
    shape = kernel.aero_lut["I"].shape

    wind_values = kernel.aero_lut["wind"].values.reshape(
        -1, 1, 1, 1, 1, 1, 1
    )

    intensities = np.broadcast_to(
        wind_values,
        shape,
    ).copy()

    kernel.aero_lut["I"] = (dims, intensities)

    return kernel

@pytest.fixture
def kernel_lut_decimal_angles(kernel_lut):
    """Create a raster with non-integer viewing angles."""

    kernel = kernel_lut

    kernel.prod.raster["sza"][:] = [
        [30.14, 30.24],
        [40.16, 40.24],
    ]

    kernel.prod.raster["vza"][:] = [
        [10.14, 10.24],
        [20.16, 20.24],
    ]

    kernel.prod.raster["raa"][:] = [
        [0.14, 0.24],
        [90.16, 90.24],
    ]

    return kernel

# ============================================================
# Tests for angles preparation
# ============================================================

class TestLutPreparationAngles:
    
    def test_truncates_angles(self, kernel_lut, mock_auxdata):
        """Check that the angle rasters are rounded as expected."""

        kernel_lut.lut_preparation()
        np.testing.assert_array_equal(
            kernel_lut.prod.raster["sza_trunc"].values, 
            [[30., 30.],[40., 40.]]
        )
        np.testing.assert_array_equal(
            kernel_lut.prod.raster["vza_trunc"].values,
            [[10., 10.], [20., 20.]],
        )

        np.testing.assert_array_equal(
            kernel_lut.prod.raster["raa_trunc"].values,
            [[0., 0.], [90., 90.]],
        )

    def test_ignores_nan_angles(self, kernel_lut, mock_auxdata):
        """Check that NaN angles are excluded from the LUT coordinates."""

        kernel_lut.prod.raster["sza"].values[0, 0] = np.nan
        kernel_lut.prod.raster["vza"].values[1, 1] = np.nan
        kernel_lut.prod.raster["raa"].values[0, 1] = np.nan

        kernel_lut.lut_preparation()

        assert not np.isnan(kernel_lut.szas).any()
        assert not np.isnan(kernel_lut.vzas).any()
        assert not np.isnan(kernel_lut.azis).any()

    def test_azimuth_values(self, kernel_lut, mock_auxdata):
        """Check the conversion from relative azimuth to LUT azimuth."""

        kernel_lut.lut_preparation()

        expected_azimuths = np.array([90., 180.])

        np.testing.assert_allclose(
            kernel_lut.azis,
            expected_azimuths,
        )
    
    def test_rounds_angles(self, kernel_lut_decimal_angles, mock_auxdata):
        """Check that angle resolution controls the rounding."""

        kernel_lut_decimal_angles.lut_preparation(
            ang_resol={
                "sza": 1,
                "vza": 1,
                "raa": 0,
            },
        )

        raster = kernel_lut_decimal_angles.prod.raster

        np.testing.assert_allclose(
            raster["sza_trunc"].values,
            [
                [30.1, 30.2],
                [40.2, 40.2],
            ],
        )

        np.testing.assert_allclose(
            raster["vza_trunc"].values,
            [
                [10.1, 10.2],
                [20.2, 20.2],
            ],
        )

        np.testing.assert_allclose(
            raster["raa_trunc"].values,
            [
                [0., 0.],
                [90., 90.],
            ],
        )

# ============================================================
# Tests for Rayleigh reflectance
# ============================================================

class TestLutPreparationRayleigh:
    
    def test_rayleigh_spectral_interpolation(self,kernel_lut_spectral,mock_auxdata):
        """Check spectral interpolation of Rayleigh reflectance."""

        kernel_lut_spectral.lut_preparation()

        result = kernel_lut_spectral.Rray

        expected_intensity = result["wl"] / 1000

        expected = (
            expected_intensity / np.cos(np.radians(result["sza"]))
        ).broadcast_like(result)

        np.testing.assert_allclose(
            result.values,
            expected.values,
            rtol=1e-6,
            atol=1e-8,
        )

    def test_rayleigh_reflectance(self, kernel_lut, mock_auxdata):
        """Check Rayleigh reflectance against explicitly expected angles."""

        kernel_lut.lut_preparation()

        result = kernel_lut.Rray

        # Rayleigh intensity is constant: I0 = 0.2.
        # The raster contains SZA values of 30 and 40 degrees.
        expected_sza = np.array([30., 40.])

        expected = 0.2 / np.cos(np.radians(expected_sza))

        np.testing.assert_allclose(
            result["sza"].values,
            expected_sza,
        )

        # The other dimensions must not affect the Rayleigh intensity.
        expected = xr.DataArray(
            expected,
            dims=("sza",),
            coords={"sza": expected_sza},
        ).broadcast_like(result)

        np.testing.assert_allclose(
            result.values,
            expected.values,
            rtol=1e-6,
        )


# ============================================================
# Tests for Aerosol model combinations
# ============================================================

class TestLutPreparationAerosol:

    def test_diffuse_reflectance_combination(self,kernel_lut_diffuse,mock_auxdata):
        """Check the weighted combination of aerosol models in Rdiff_lut."""

        weights = np.array([0., 0.5, 0., 0.5, 0., 0.])

        kernel_lut_diffuse.lut_preparation(
            aerosol_combination=weights,
        )

        result = kernel_lut_diffuse.Rdiff_lut

        # Models have constant intensities: 1, 2, 3, 4, 5, 6.
        # The selected weights give: 0.5 * 2 + 0.5 * 4 = 3.
        expected = 3. / np.cos(np.radians(result["sza"]))
        expected = expected.broadcast_like(result)

        assert len(kernel_module.AEROSOL_MODELS) == 6
        
        np.testing.assert_allclose(
            result.values,
            expected.values,
            rtol=1e-6,
        )

    def test_invalid_aerosol_combination(self,kernel_lut,mock_auxdata,caplog):
        """Check that invalid aerosol weights stop LUT preparation."""

        with caplog.at_level(logging.ERROR):
            result = kernel_lut.lut_preparation(
                aerosol_combination="invalid",
            )

        assert result is None

        assert any(
            "aerosol_combination" in record.message
            for record in caplog.records
        )

        assert not hasattr(kernel_lut, "Rdiff_lut")

    def test_accepts_dataarray_combination(self,kernel_lut,mock_auxdata):
        """Check that aerosol weights can be supplied as a DataArray."""

        weights = xr.DataArray(
            [0., 0.5, 0., 0.5, 0., 0.],
            dims=("model",),
            coords={"model": kernel_module.AEROSOL_MODELS},
        )

        kernel_lut.lut_preparation(
            aerosol_combination=weights,
        )

        xr.testing.assert_equal(
            kernel_lut.aerosol_combination,
            weights,
        )

    def test_default_aerosol_combination(self,kernel_lut,mock_auxdata):
        """Check that the default aerosol combination is used."""

        kernel_lut.lut_preparation()

        xr.testing.assert_equal(
            kernel_lut.aerosol_combination,
            xr.DataArray(
                kernel_module.AEROSOL_COMBINATION,
                coords={"model": kernel_module.AEROSOL_MODELS},
            ),
        )

# ============================================================
# Tests for AOT preparation
# ============================================================

class TestLutPreparationAOT:

    def test_combination_and_interpolation(self,kernel_lut_aot,mock_auxdata):
        """Check weighted AOT combination and interpolation."""

        weights = np.array([0., 0.5, 0., 0.5, 0., 0.])

        kernel_lut_aot.lut_preparation(
            aerosol_combination=weights,
        )

        result = kernel_lut_aot.aot_lut

        expected = 3.0 * (
            result["aot_ref"] + result["wl"] / 1000
        )

        expected = expected.broadcast_like(result)

        assert len(kernel_module.AEROSOL_MODELS) == 6
        np.testing.assert_allclose(
            result.values,
            expected.values,
            rtol=1e-6,
            atol=1e-8,
        )
    
    def test_aot_reference_coordinates(self,kernel_lut,mock_auxdata):
        """Check the reference AOT coordinates after interpolation."""

        requested_aot_refs = np.linspace(0.0, 0.8, 25)

        kernel_lut.lut_preparation(
            aot_refs=requested_aot_refs,
        )

        np.testing.assert_allclose(
            kernel_lut.aot_refs,
            requested_aot_refs,
        )

        np.testing.assert_allclose(
            kernel_lut.Rdiff_lut["aot_ref"].values,
            requested_aot_refs,
        )

        np.testing.assert_allclose(
            kernel_lut.aot_lut["aot_ref"].values,
            requested_aot_refs,
        )

        np.testing.assert_allclose(
            kernel_lut.trans_aero_lut["aot_ref"].values,
            requested_aot_refs,
        )


# ============================================================
# Tests for Transmittance preparation
# ============================================================

class TestLutPreparationTransmittance:

    def test_transmittance_combination(self,kernel_lut_transmittance,mock_auxdata):
        """Check the weighted combination of transmittance models."""

        weights = np.array([0., 0.5, 0., 0.5, 0., 0.])

        kernel_lut_transmittance.lut_preparation(
            aerosol_combination=weights,
        )

        result = kernel_lut_transmittance.trans_aero_lut["Ttot_Ed"]

        # The selected models have transmittances 2 and 4.
        # Their weighted mean is 3.
        expected = xr.full_like(result, 3.)

        np.testing.assert_allclose(
            result.values,
            expected.values,
            rtol=1e-6,
            atol=1e-8,
        )

# ============================================================
# Tests for LUT ouput and structure
# ============================================================

class TestLutPreparationOutputs:

    def test_removes_model_dimension(self,kernel_lut,mock_auxdata):
        """Check that no prepared LUT retains the aerosol model dimension."""

        kernel_lut.lut_preparation()

        for name in (
            "Rray",
            "Rdiff_lut",
            "aot_lut",
            "trans_aero_lut",
        ):
            result = getattr(kernel_lut, name)

            assert "model" not in result.dims, (
                f"{name} unexpectedly contains the model dimension"
            )

    def test_loads_prepared_luts(self,kernel_lut,mock_auxdata,monkeypatch):
        """Check that _load_small is called for all prepared LUTs."""

        original = kernel_module._load_small
        calls = []

        def spy_load_small(obj, *args, **kwargs):
            calls.append(obj)
            return original(obj, *args, **kwargs)

        monkeypatch.setattr(
            kernel_module,
            "_load_small",
            spy_load_small,
        )

        kernel_lut.lut_preparation()

        # Four LUTs and two auxiliary datasets.
        assert len(calls) >= 4

        assert kernel_lut.trans_aero_lut is not None
        assert kernel_lut.Rray is not None
        assert kernel_lut.Rdiff_lut is not None
        assert kernel_lut.aot_lut is not None

    def test_output_dimensions(self,kernel_lut, mock_auxdata):
        """Check dimensions and wavelength coordinates of prepared LUTs."""

        kernel_lut.lut_preparation()

        expected_wavelengths = kernel_lut.wl_true

        for result in (
            kernel_lut.Rray,
            kernel_lut.Rdiff_lut,
            kernel_lut.aot_lut,
            kernel_lut.trans_aero_lut,
        ):
            np.testing.assert_allclose(
                result["wl"].values,
                expected_wavelengths,
            )

        assert "model" not in kernel_lut.Rdiff_lut.dims
        assert "model" not in kernel_lut.aot_lut.dims
        assert "model" not in kernel_lut.trans_aero_lut.dims

        assert "aot_ref" in kernel_lut.Rdiff_lut.dims
        assert "aot_ref" in kernel_lut.aot_lut.dims
        assert "aot_ref" in kernel_lut.trans_aero_lut.dims

        np.testing.assert_allclose(
            kernel_lut.aot_refs,
            kernel_lut.Rdiff_lut["aot_ref"].values,
        )

        np.testing.assert_allclose(
            kernel_lut.szas,
            kernel_lut.Rdiff_lut["sza"].values,
        )

        np.testing.assert_allclose(
            kernel_lut.vzas,
            kernel_lut.Rdiff_lut["vza"].values,
        )

        np.testing.assert_allclose(
            kernel_lut.azis,
            kernel_lut.Rdiff_lut["azi"].values,
        )

# ============================================================
# Tests for preparation caching
# ============================================================

class TestLutPreparationCaching:
        
    def test_uses_separate_cache_for_angle_resolution(self,kernel_lut,mock_auxdata,monkeypatch):
        """Check that different angle resolutions use different cache entries."""

        original = kernel_module._unique_values
        calls = []

        def spy_unique_values(*arrays):
            calls.append(arrays)
            return original(*arrays)

        monkeypatch.setattr(
            kernel_module,
            "_unique_values",
            spy_unique_values,
        )

        kernel_lut.lut_preparation()

        kernel_lut.lut_preparation(
            ang_resol={
                "sza": 0,
                "vza": 0,
                "raa": 0,
            },
        )

        assert len(calls) == 2

    def test_caches_unique_angles(self,kernel_lut,mock_auxdata,monkeypatch):
        """Check that unique angle values are cached."""

        original = kernel_module._unique_values
        calls = []

        def spy_unique_values(*arrays):
            calls.append(arrays)
            return original(*arrays)

        monkeypatch.setattr(
            kernel_module,
            "_unique_values",
            spy_unique_values,
        )

        kernel_lut.lut_preparation()

        assert len(calls) == 1

        # Repeating the preparation with the same angular resolution
        # must reuse the cached angle values.
        kernel_lut.lut_preparation()

        assert len(calls) == 1

# ============================================================
# Tests for preparation auxiliary data and configuration
# ============================================================

class TestLutPreparationAuxiliaryData:

    def test_selects_nearest_wind(self,kernel_lut_wind_dependent, mock_auxdata):
        """Check that the nearest wind-dependent LUT is selected."""

        kernel_lut_wind_dependent.lut_preparation(
            wind=2.4,
        )

        result = kernel_lut_wind_dependent.Rdiff_lut

        expected = 2. / np.cos(np.radians(result["sza"]))
        expected = expected.broadcast_like(result)

        np.testing.assert_allclose(
            result.values,
            expected.values,
            rtol=1e-6,
        )

    def test_auxiliary_data(self,kernel_lut, mock_auxdata):
        """Check that auxiliary spectral data are initialized."""

        kernel_lut.lut_preparation()

        np.testing.assert_allclose(
            kernel_lut.sunglint_eps.values,
            np.ones(len(kernel_lut.wl_true)),
        )

        np.testing.assert_allclose(
            kernel_lut.rot.values,
            np.full(len(kernel_lut.wl_true), 0.1),
        )

        np.testing.assert_allclose(
            kernel_lut.sunglint_eps["wl"].values,
            kernel_lut.wl_true,
        )

        np.testing.assert_allclose(
            kernel_lut.rot["wl"].values,
            kernel_lut.wl_true,
        )








