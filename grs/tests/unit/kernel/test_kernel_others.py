from unittest.mock import PropertyMock, patch

import numpy as np
import pytest
import xarray as xr

from types import SimpleNamespace

from grs.grs_kernel import Kernel

# ============================================================
# Fixtures
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
# Tests mu_N
# ============================================================

class TestMu_N:
     
    def test_monoview(self):
        """Check mu_N for zero angles in monoview mode."""
        kernel = Kernel.__new__(Kernel)
        sza = xr.DataArray(
                np.full((2, 2), 60.),
                dims=("y", "x"),
        )

        vza = xr.DataArray(
            np.full((2, 2), 60.),
            dims=("y", "x"),
        )

        azi = xr.DataArray(
            np.full((2, 2), 180.),
            dims=("y", "x"),
        )

        result = kernel.mu_N(
            sza=sza,
            vza=vza,
            azi=azi,
            monoview=True,
        )

        np.testing.assert_allclose(
            result.values,
            np.ones((2, 2)),
            rtol=1e-12,
            atol=1e-12,
        )
        
    def test_multiview(self,kernel):
        """Check mu_N calculation and dimensions in multiview mode."""

        wl = [500., 600., 700., 800.]

        # Angles vary with wavelength.
        sza_values = np.array([0., 30., 45., 60.])
        vza_values = np.array([0., 10., 20., 30.])

        shape = (4, 2, 2)

        sza = xr.DataArray(
            np.broadcast_to(sza_values[:, None, None], shape).copy(),
            dims=("wl", "y", "x"),
            coords={"wl": wl},
        )
        
        vza = xr.DataArray(
            np.broadcast_to(vza_values[:, None, None], shape).copy(),
            dims=("wl", "y", "x"),
            coords={"wl": wl},
        )

        azi = xr.DataArray(
            np.zeros(shape),
            dims=("wl", "y", "x"),
            coords={"wl": wl},
        )

        result = kernel.mu_N(
            sza=sza,
            vza=vza,
            azi=azi,
            monoview=False,
        )

        assert result.dims == ("wl", "y", "x")

        expected_values = np.cos(
            (np.radians(sza_values) + np.radians(vza_values)) / 2
        )

        expected = np.broadcast_to(
            expected_values[:, None, None],
            shape,
        )

        np.testing.assert_allclose(
            result.values,
            expected,
            rtol=1e-12,
            atol=1e-12,
        )


# ============================================================
# Tests conv_mapping
# ============================================================

class TestConvMapping:
     
    def test_ignores_nan_when_computing_mean(self):
        kernel = Kernel.__new__(Kernel)
        test_array = np.array([2.0,np.nan,4.0])
        result = kernel.conv_mapping(test_array)
        
        assert result == 3.0

    def test_returns_central_value_when_neighbours_are_nan(self):
        kernel = Kernel.__new__(Kernel)
        test_array = np.array([np.nan, 3.0, np.nan])
        result = kernel.conv_mapping(test_array)
        
        assert result == 3.0

    def test_finite_values(self):
        kernel = Kernel.__new__(Kernel)
        test_array = np.array([2.0, 3.0, 4.0])
        result = kernel.conv_mapping(test_array)
        
        assert result == 3.0


# ============================================================
# Tests smoothing
# ============================================================

class TestSmoothing:
     
    def test_uses_nanmean_box(self,kernel, monkeypatch):
        raster = xr.Dataset(
            {
                "wind": (
                    ("y", "x"),
                    np.array([
                        [1.0, 2.0],
                        [3.0, 4.0],
                    ]),
                )
            },
            coords={
                "y": [10, 20],
                "x": [100, 200],
            },
        )

        filtered = np.array([
            [1.5, 2.5],
            [3.5, 4.5],
        ])

        expected = np.array([
            [2.0, 2.5],
            [3.5, 4.0],
        ])

        calls = {}

        def fake_filter2d(image, weights, windows):
            calls["filter2d"] = True
            np.testing.assert_array_equal(image, raster["wind"].values)
            np.testing.assert_array_equal(weights, np.ones((2, 2)))
            np.testing.assert_array_equal(windows, [3, 3])
            return filtered

        def fake_nanmean_box(image, size):
            calls["nanmean_box"] = True
            np.testing.assert_array_equal(image, filtered)
            assert size == (3, 3)
            return expected

        monkeypatch.setattr(kernel, "filter2d", fake_filter2d)
        monkeypatch.setattr(kernel, "_nanmean_box", fake_nanmean_box)

        result = kernel.smoothing(
            raster,
            windows=np.array([3, 3]),
            mask=np.ones((3, 3)),
        )

        assert calls["filter2d"]
        assert calls["nanmean_box"]

        np.testing.assert_array_equal(result["wind"].values, expected)
        np.testing.assert_array_equal(result["y"].values, [10, 20])
        np.testing.assert_array_equal(result["x"].values, [100, 200])

    def test_uses_generic_filter(self,monkeypatch):
        kernel = Kernel.__new__(Kernel)
        raster = xr.Dataset(
            {
                "wind": (
                    ("y", "x"),
                    np.array([
                        [1.0, 2.0],
                        [3.0, 4.0],
                    ]),
                )
            },
            coords={
                "y": [10, 20],
                "x": [100, 200],
            },
        )

        filtered = np.array([
            [1.5, 2.5],
            [3.5, 4.5],
        ])

        expected = np.array([
            [2.0, 2.5],
            [3.5, 4.0],
        ])

        calls = {}

        monkeypatch.setattr(
            kernel,
            "filter2d",
            lambda image, weights, windows: filtered,
        )

        def fail_nanmean_box(*args, **kwargs):
            pytest.fail("_nanmean_box ne devrait pas être appelé")

        monkeypatch.setattr(kernel, "_nanmean_box", fail_nanmean_box)

        def fake_generic_filter(image, function, footprint, mode):
            calls["function"] = function
            calls["footprint"] = footprint
            calls["mode"] = mode
            np.testing.assert_array_equal(image, filtered)
            return expected

        monkeypatch.setattr(
            "scipy.ndimage.generic_filter",
            fake_generic_filter,
        )

        mask = np.array([
            [0, 1, 0],
            [1, 1, 1],
            [0, 1, 0],
        ])

        result = kernel.smoothing(
            raster,
            windows=np.array([3, 3]),
            mask=mask,
        )

        assert calls["function"] == Kernel.conv_mapping
        assert calls["mode"] == "nearest"
        np.testing.assert_array_equal(calls["footprint"], mask)
        np.testing.assert_array_equal(result["wind"].values, expected)

# ============================================================
# Tests _nanmean_box
# ============================================================

class TestNanMeanBox:
     
    @pytest.mark.parametrize(
        "image, size, expected",
        [
            # local mean withou NaN
            (
                np.array([[1., 2., 3.],
                        [4., 5., 6.],
                        [7., 8., 9.]]),
                (3, 3),
                np.array([
                    [7/3, 3., 11/3],
                    [13/3, 5., 17/3],
                    [19/3, 7., 23/3]
                ])
            ),
            # Horizontal window
            (
                np.array([[1., 2., 3., 4.]]),
                (1, 3),
                np.array([[4/3, 2., 3., 11/3]])
            ),
            # Size 1 window: values remain unchanged 
            (
                np.array([[1., np.nan],
                        [3., 4.]]),
                (1, 1),
                np.array([[1., np.nan],
                        [3., 4.]])
            ),
            # Array full of NaN
            (
                np.full((3, 3), np.nan),
                (3, 3),
                np.full((3, 3), np.nan)
            ),
        ],
    )
    def test_nominal(self,image, size, expected):
        result = Kernel._nanmean_box(image, size)
        np.testing.assert_allclose(result, expected, equal_nan=True)

    def test_single_valid_pixel(self):
        """ONly one valid pixel should give its value to the whole neighbourhoud"""
        image = np.full((3, 3), np.nan)
        image[1, 1] = 42.

        result = Kernel._nanmean_box(image, (3, 3))

        np.testing.assert_allclose(result, 42.)

    def test_nearest_edge_extension(self):
        """Checks the process on edges with mode='nearest'"""
        image = np.array([[1., 2.],
                        [3., 4.]])

        result = Kernel._nanmean_box(image, (3, 3))

        expected = np.array([[2., 7/3],
                            [8/3, 3.]])

        np.testing.assert_allclose(result, expected)

    def test_output_properties(self):
        """The ouptut maintains the shape and uses float64"""
        image = np.array([[1, 2], [3, 4]], dtype=np.int32)

        result = Kernel._nanmean_box(image, (3, 3))

        assert result.shape == image.shape
        assert result.dtype == np.float64


# ============================================================
# Tests snapshot
# ============================================================

class TestSnapshot:

    def test_copies_state_and_drops_heavy_attributes(self):
        """Test that _snapshot creates a lightweight copy of the kernel."""
        kernel = Kernel.__new__(Kernel)

        # Set up the product and the attributes required by _snapshot.
        kernel.prod = SimpleNamespace(_type=np.float32, sensor='S2A')
        kernel._cache = {'cached_value': np.array([1, 2, 3])}

        kernel.wind_img = xr.Dataset({
            'wind': (('y', 'x'), np.array([[2.0, 3.0], [4.0, 5.0]]))
        })

        # Simulate the heavy attributes that must not be retained.
        kernel.raster = object()
        kernel.aero_lut = object()
        kernel.trans_lut = object()
        kernel.cams = object()
        kernel.coarse_masked_raster = object()
        kernel.xres = object()
        kernel.aot_ref_img = object()
        kernel.aot_ref_raster = object()

        # Mock the property to avoid loading actual surface pressure data.
        surface_pressure = xr.DataArray(
            [[101325.0, 101300.0], [101200.0, 101100.0]],
            dims=('y', 'x'),
        )

        with patch.object(
            Kernel,
            '_surface_pressure',
            new_callable=PropertyMock,
            create=True,
            return_value=surface_pressure,
        ) as mocked_pressure:
            state = kernel._snapshot()

        # Verify that the original kernel is not modified.
        for name in (
            'raster',
            'aero_lut',
            'trans_lut',
            'cams',
            'coarse_masked_raster',
            'xres',
            'aot_ref_img',
            'aot_ref_raster',
        ):
            assert hasattr(kernel, name)

        # Verify that heavy attributes are removed from the snapshot.
        for name in (
            'raster',
            'aero_lut',
            'trans_lut',
            'cams',
            'coarse_masked_raster',
            'xres',
            'aot_ref_img',
            'aot_ref_raster',
        ):
            assert not hasattr(state, name)

        # Verify that lightweight attributes are preserved.
        assert state is not kernel
        assert state.prod._type is np.float32
        assert state.prod.sensor == 'S2A'
        assert state._cache == kernel._cache
        assert state._cache is not kernel._cache

        # Verify that the wind dataset is copied without deep-copying its data.
        assert state.wind_img is not kernel.wind_img
        xr.testing.assert_identical(state.wind_img, kernel.wind_img)

        # Verify that surface pressure is accessed before the CAMS reference is removed.
        mocked_pressure.assert_called_once_with()