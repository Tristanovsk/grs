"""Unit tests kernel.set_gas_transmittance"""

import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel
import grs.grs_kernel as kernel_module

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def kernel_for_gas_transmittance():
    """Create a Kernel instance without running its constructor."""
    kernel = Kernel.__new__(Kernel)
    kernel.prod = object()
    kernel.cams = object()
    return kernel


@pytest.fixture
def gas_transmittance_data():
    """Create synthetic gaseous transmittance rasters."""
    coords = {
        "wl": [443, 560, 665],
        "y": [10, 20],
        "x": [100, 200],
    }

    diffuse = xr.DataArray(
        np.full((3, 2, 2), 0.8, dtype=np.float32),
        dims=("wl", "y", "x"),
        coords=coords,
        name="transmittance",
    )

    total = xr.DataArray(
        np.full((3, 2, 2), 0.9, dtype=np.float32),
        dims=("wl", "y", "x"),
        coords=coords,
        name="transmittance",
    )

    return diffuse, total

# ============================================================
# Tests Gas Transmittance initialization
# ============================================================

class TestSetGasTransmittanceInitialization:
    
    def test_initializes_both_rasters(self,monkeypatch,kernel_for_gas_transmittance,gas_transmittance_data):
        """Test that diffuse and total gaseous transmittance rasters are initialized."""
        diffuse, total = gas_transmittance_data
        calls = []

        class FakeGaseousTransmittance:
            def __init__(self, prod, cams):
                assert prod is kernel_for_gas_transmittance.prod
                assert cams is kernel_for_gas_transmittance.cams

                self.coef_abs_scat = {
                    "h2o": 1,
                    "ch4": 1,
                }

            def get_gaseous_transmittance(self, gases, background):
                calls.append({
                    "gases": list(gases),
                    "background": background,
                    "coefficients": dict(self.coef_abs_scat),
                })

                if self.coef_abs_scat["h2o"] == 0.5:
                    return diffuse.transpose("x", "wl", "y")

                return total.transpose("y", "x", "wl")

        loaded = []

        def fake_load_small(data):
            loaded.append(data)
            return data

        monkeypatch.setattr(
            kernel_module,
            "GaseousTransmittance",
            FakeGaseousTransmittance,
        )
        monkeypatch.setattr(kernel_module, "_load_small", fake_load_small)

        kernel = kernel_for_gas_transmittance
        kernel.set_gas_transmittance()

        assert hasattr(kernel, "Tg_diff_raster")
        assert hasattr(kernel, "Tg_raster")

        xr.testing.assert_equal(kernel.Tg_diff_raster, diffuse)
        xr.testing.assert_equal(kernel.Tg_raster, total)

        assert len(loaded) == 2
        assert len(calls) == 2

    def test_calls_load_small_twice(self,monkeypatch,kernel_for_gas_transmittance,gas_transmittance_data):
        """Test that both transmittance rasters are passed through the loading helper."""
        diffuse, total = gas_transmittance_data
        loaded = []

        class FakeGaseousTransmittance:
            def __init__(self, prod, cams):
                self.coef_abs_scat = {
                    "h2o": 1,
                    "ch4": 1,
                }

            def get_gaseous_transmittance(self, gases, background):
                if self.coef_abs_scat["h2o"] == 0.5:
                    return diffuse

                return total

        def fake_load_small(data):
            loaded.append(data.copy())
            return data

        monkeypatch.setattr(
            kernel_module,
            "GaseousTransmittance",
            FakeGaseousTransmittance,
        )
        monkeypatch.setattr(kernel_module, "_load_small", fake_load_small)

        kernel_for_gas_transmittance.set_gas_transmittance()

        assert len(loaded) == 2
        assert loaded[0].dims == ("wl", "y", "x")
        assert loaded[1].dims == ("wl", "y", "x")

# ============================================================
# Tests Gas Transmittance configuration
# ============================================================

class TestSetGasTransmittanceConfiguration:
    
    def test_uses_expected_gases_and_coefficients(self,monkeypatch,kernel_for_gas_transmittance,gas_transmittance_data):
        """Test the gases and absorption-scattering coefficients used in each calculation."""
        diffuse, total = gas_transmittance_data
        calls = []

        class FakeGaseousTransmittance:
            def __init__(self, prod, cams):
                self.coef_abs_scat = {
                    "h2o": 1,
                    "ch4": 1,
                }

            def get_gaseous_transmittance(self, gases, background):
                calls.append({
                    "gases": list(gases),
                    "background": background,
                    "coefficients": dict(self.coef_abs_scat),
                })

                if len(calls) == 1:
                    return diffuse

                return total

        monkeypatch.setattr(
            kernel_module,
            "GaseousTransmittance",
            FakeGaseousTransmittance,
        )
        monkeypatch.setattr(kernel_module, "_load_small", lambda data: data)

        kernel_for_gas_transmittance.set_gas_transmittance()

        assert calls[0]["gases"] == ["h2o", "ch4"]
        assert calls[0]["background"] is False
        assert calls[0]["coefficients"]["h2o"] == 0.5
        assert calls[0]["coefficients"]["ch4"] == 0.65

        assert calls[1]["gases"] == ["h2o", "ch4"]
        assert calls[1]["background"] is False
        assert calls[1]["coefficients"]["h2o"] == 1
        assert calls[1]["coefficients"]["ch4"] == 1

# ============================================================
# Tests Gas Transmittance outputs
# ============================================================

class TestSetGasTransmittanceOutputs:
    
    def test_transposes_dimensions(self,monkeypatch,kernel_for_gas_transmittance,gas_transmittance_data):
        """Test that both output rasters use the (wl, y, x) dimension order."""
        diffuse, total = gas_transmittance_data

        class FakeGaseousTransmittance:
            def __init__(self, prod, cams):
                self.coef_abs_scat = {
                    "h2o": 1,
                    "ch4": 1,
                }

            def get_gaseous_transmittance(self, gases, background):
                if self.coef_abs_scat["h2o"] == 0.5:
                    return diffuse.transpose("x", "y", "wl")

                return total.transpose("y", "x", "wl")

        monkeypatch.setattr(
            kernel_module,
            "GaseousTransmittance",
            FakeGaseousTransmittance,
        )
        monkeypatch.setattr(kernel_module, "_load_small", lambda data: data)

        kernel_for_gas_transmittance.set_gas_transmittance()

        assert kernel_for_gas_transmittance.Tg_diff_raster.dims == ("wl", "y", "x")
        assert kernel_for_gas_transmittance.Tg_raster.dims == ("wl", "y", "x")

    def test_preserves_transmittance_values(self,monkeypatch,kernel_for_gas_transmittance,gas_transmittance_data):
        """Test that the loading helper does not alter transmittance values."""
        diffuse, total = gas_transmittance_data

        class FakeGaseousTransmittance:
            def __init__(self, prod, cams):
                self.coef_abs_scat = {
                    "h2o": 1,
                    "ch4": 1,
                }

            def get_gaseous_transmittance(self, gases, background):
                if self.coef_abs_scat["h2o"] == 0.5:
                    return diffuse

                return total

        monkeypatch.setattr(
            kernel_module,
            "GaseousTransmittance",
            FakeGaseousTransmittance,
        )
        monkeypatch.setattr(kernel_module, "_load_small", lambda data: data)

        kernel_for_gas_transmittance.set_gas_transmittance()

        np.testing.assert_allclose(
            kernel_for_gas_transmittance.Tg_diff_raster.values,
            0.8,
        )
        np.testing.assert_allclose(
            kernel_for_gas_transmittance.Tg_raster.values,
            0.9,
        )