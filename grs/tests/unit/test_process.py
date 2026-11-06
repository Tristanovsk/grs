import os
from unittest.mock import Mock

import numpy as np
import pytest
import xarray as xr
import dask.array as da

import grs.grs_process as grs_process

from grs.grs_process import _correct_gas_absorption, _interp_like, _cams_from_archive


def test_correct_gas_absorption_no_chunks():
    x_bands = np.array([0., 1., 2.])
    y_bands = np.array([0., 1., 2.])

    wl = np.array([500., 600.])

    bands_values = np.array([
        # Bande 500 nm
        [
            [0.4, 0.5, 0.6],
            [0.5, 0.6, 0.7],
            [0.6, 0.7, 0.8],
        ],
        # Bande 600 nm
        [
            [0.2, 0.3, 0.4],
            [0.3, 0.4, 0.5],
            [0.4, 0.5, 0.6],
        ],
    ])

    bands = xr.DataArray(
        bands_values.copy(),
        dims=("wl", "y", "x"),
        coords={
            "wl": wl,
            "y": y_bands,
            "x": x_bands,
        },
    )
    
    x_tg = np.array([0., 2.])
    y_tg = np.array([0., 2.])
    
    Tg_values = np.array([
        [
            [0.5, 0.7],
            [0.9, 1.1],
        ],
        [
            [0.8, 0.7],
            [0.6, 0.5],
        ],
    ])

    Tg_raster = xr.DataArray(
        Tg_values,
        dims=("wl", "y", "x"),
        coords={
            "wl": wl,
            "y": y_tg,
            "x": x_tg,
        },
    )
    
    Tg_expected = np.array([
        [
            [0.5, 0.6, 0.7],
            [0.7, 0.8, 0.9],
            [0.9, 1.0, 1.1],
        ],
        [
            [0.8, 0.75, 0.7],
            [0.7, 0.65, 0.6],
            [0.6, 0.55, 0.5],
        ],
    ])
    
    expected = bands_values / Tg_expected
    
    result = _correct_gas_absorption(bands=bands, Tg_raster=Tg_raster)
    xr.testing.assert_allclose(
        result,
        xr.DataArray(
            expected,
            dims=("wl", "y", "x"),
            coords={
                "wl": wl,
                "y": y_bands,
                "x": x_bands,
            },
        ),
        rtol=1e-7,
        atol=1e-7,
    )

# ============================================================
# _interp_like
# ============================================================

class TestInterpLike:
    
    def test_nominal(self):
        coarse = xr.DataArray(
            np.array([
                [0, 20],
                [40, 60],
            ]),
            dims=("y", "x"),
            coords={
                "y": np.array([0, 2]),
                "x": np.array([0, 2]),
                "time": np.datetime64("2025-01-01"),
            },
            name= "original"
        )
        
        like = xr.DataArray(
            da.zeros((5, 5), chunks=(2, 2)),
            dims=("y", "x"),
            coords={
                "y": np.arange(5) / 2,
                "x": np.arange(5) / 2,
                "quality": (
                    ("y","x"),
                    np.ones((5,5)),
                )
            },
        )
        
        result = _interp_like(coarse, like, name="interpolated")
        computed = result.compute()
        expected = np.array([
            [0, 5, 10, 15,  20],
            [10, 15, 20, 25, 30],
            [20, 25, 30, 35, 40],
            [30, 35, 40, 45, 50],
            [40, 45, 50, 55, 60],
        ])
        assert isinstance(result.data, da.Array)
        assert result.shape == like.shape
        assert result.dims == ("y", "x")
        assert result.dtype == np.float64
        assert "x" in result.coords
        assert "y" in result.coords
        assert "time" in result.coords
        assert result.coords["time"].ndim == 0
        np.testing.assert_equal(
            result.coords["time"].values,
            np.datetime64("2025-01-01"),
        )
        assert result.name == "interpolated"
        assert "quality" not in result.coords
        np.testing.assert_allclose(result, expected)

    def test_outside_source_bounds(self):
        coarse = xr.DataArray(
            [[0.0, 1.0], [2.0, 3.0]],
            dims=("y", "x"),
            coords={
                "y": [0, 2],
                "x": [0, 2],
            },
        )

        like = xr.DataArray(
            np.zeros((3, 4)),
            dims=("y", "x"),
            coords={
                "y": [0, 1, 2],
                "x": [0, 1, 2, 3],
            },
        )

        result = _interp_like(coarse, like)

        assert result.sel(x=3).isnull().all()

# ============================================================
# _cams_from_archive
# ============================================================

class TestCamsFromArchive:
    
    def test_file_exists(self,monkeypatch, caplog):
        raster = Mock()
        raster.attrs = {"tile": "T31ABC"}

        cams = Mock()
        cams.filepath = "/cams/2025/01/01/cams.nc"
        cams.date_str = "2025-01-01"

        constructor_calls = []

        def fake_cams_product(*args, **kwargs):
            constructor_calls.append((args, kwargs))
            return cams

        monkeypatch.setattr(grs_process, "CamsProduct", fake_cams_product)

        monkeypatch.setattr(
            grs_process.os.path,
            "exists",
            lambda path: path == cams.filepath,
        )

        result = _cams_from_archive(raster, "/cams")

        assert result is cams

        assert len(constructor_calls) == 1
        assert constructor_calls[0][0] == (raster,)
        assert constructor_calls[0][1] == {
            "dir": "/cams",
            "suffix": "_T31ABC",
        }

    def test_uses_legacy_file(self,monkeypatch):
        raster = Mock()
        raster.attrs = {"tile": "T31ABC"}

        cams_archive = Mock()
        cams_archive.filepath = "/cams/archive.nc"
        cams_archive.date_str = "2025-01-01"

        cams_legacy = Mock()
        cams_legacy.filepath = "/cams/T31ABC/2025-01-01-forecast_T31ABC.nc"

        constructor_calls = []

        def fake_cams_product(*args, **kwargs):
            constructor_calls.append((args, kwargs))

            if "cams_file" in kwargs:
                return cams_legacy

            return cams_archive

        monkeypatch.setattr(grs_process, "CamsProduct", fake_cams_product)

        legacy_path = os.path.join(
            "/cams",
            "T31ABC",
            "2025-01-01-forecast_T31ABC.nc",
        )

        existing_files = {legacy_path}

        monkeypatch.setattr(
            grs_process.os.path,
            "exists",
            lambda path: path in existing_files,
        )

        result = _cams_from_archive(raster, "/cams")

        assert result is cams_legacy

        assert len(constructor_calls) == 2

        assert constructor_calls[1][0] == (raster,)
        assert constructor_calls[1][1] == {
            "cams_file": legacy_path,
        }

    def test_no_existing_file(self,monkeypatch):
        raster = Mock()
        raster.attrs = {"tile": "T31ABC"}

        cams = Mock()
        cams.filepath = "/cams/archive.nc"
        cams.date_str = "2025-01-01"

        constructor_calls = []

        def fake_cams_product(*args, **kwargs):
            constructor_calls.append((args, kwargs))
            return cams

        monkeypatch.setattr(grs_process, "CamsProduct", fake_cams_product)

        monkeypatch.setattr(
            grs_process.os.path,
            "exists",
            lambda path: False,
        )

        result = _cams_from_archive(raster, "/cams")

        assert result is cams
        assert len(constructor_calls) == 1

# ============================================================
# _aerosol_settings
# ============================================================

class TestAerosolSettings:
    
    def test_from_config(self,monkeypatch):
        config = {
            key: f"value_{key}"
            for key in grs_process.AEROSOL_SETTINGS
        }

        monkeypatch.setattr(
            grs_process,
            "load_config",
            lambda: {"aerosol": config},
        )

        result = grs_process._aerosol_settings()

        assert result == config

    def test_override(self,monkeypatch):
        config = {
            key: f"value_{key}"
            for key in grs_process.AEROSOL_SETTINGS
        }

        key = grs_process.AEROSOL_SETTINGS[0]

        monkeypatch.setattr(
            grs_process,
            "load_config",
            lambda: {"aerosol": config},
        )

        result = grs_process._aerosol_settings(**{key: "new_value"})

        assert result[key] == "new_value"

    def test_none_override(self,monkeypatch):
        config = {
            key: f"value_{key}"
            for key in grs_process.AEROSOL_SETTINGS
        }

        key = grs_process.AEROSOL_SETTINGS[0]

        monkeypatch.setattr(
            grs_process,
            "load_config",
            lambda: {"aerosol": config},
        )

        result = grs_process._aerosol_settings(**{key: None})

        assert result[key] == config[key]

    def test_missing_key(self,monkeypatch):
        missing_key = grs_process.AEROSOL_SETTINGS[0]

        config = {
            key: f"value_{key}"
            for key in grs_process.AEROSOL_SETTINGS
            if key != missing_key
        }

        monkeypatch.setattr(
            grs_process,
            "load_config",
            lambda: {"aerosol": config},
        )

        with pytest.raises(KeyError, match="aerosol settings missing"):
            grs_process._aerosol_settings()
        
    def test_missing_key_provided_by_override(self,monkeypatch):
        missing_key = grs_process.AEROSOL_SETTINGS[0]

        config = {
            key: f"value_{key}"
            for key in grs_process.AEROSOL_SETTINGS
            if key != missing_key
        }

        monkeypatch.setattr(
            grs_process,
            "load_config",
            lambda: {"aerosol": config},
        )

        result = grs_process._aerosol_settings(**{missing_key: "provided_value"})

        assert result[missing_key] == "provided_value"
