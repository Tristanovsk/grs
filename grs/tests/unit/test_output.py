import builtins
import logging
import shutil
import sys
import types

import numpy as np 
import pytest 
import xarray as xr

import grs.output as mod
from unittest.mock import MagicMock, patch
from grs.output import ( L2aProduct, _base_encoding, _clear_store, _coarsen, _compression_encoding, _consolidate, _pyramid_levels, _rrs_encoding, _coord_encoding, _package_version, _reset_encoding, _load_coords, _format_from_name, _is_remote, _update_geotransform, _var_chunking, _xarray_supports_shards, _zarr_major)
from grs.output import (NETCDF_ENGINE, _CONFLICTING_ATTRS, _NETCDF_EXTENSIONS, _ZARR_EXTENSIONS, _FORMAT_EXTENSION, PYRAMID_FACTOR, PYRAMID_MIN_SIZE, LEVEL_GROUP, ANCILLARY_GROUP, _DEFAULT_SHARDS, _DEFAULT_CHUNKS_SHARDED, _DEFAULT_CHUNKS, )

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture 
def sample_dataset(): 
    """Crée un Dataset avec des encodages préexistants."""

    ds = xr.Dataset(
        data_vars={
            "Rrs": (
                ("y", "x"),
                np.array(
                    [[0.1, 0.2], [0.3, 0.4]],
                    dtype=np.float32,
                ),
            ),
            "flags": (
                ("y", "x"),
                np.array(
                    [[0, 1], [2, 3]],
                    dtype=np.uint8,
                ),
            ),
        },
        coords={
            "x": np.array([10.0, 11.0]),
            "y": np.array([45.0, 44.0]),
            "time": np.datetime64("2024-06-15"),
        },
    )

    # Simule les encodages hérités d'un fichier d'entrée
    ds["Rrs"].encoding = {
        "dtype": "int16",
        "scale_factor": 0.0001,
        "_FillValue": -9999,
        "chunksizes": (128, 128),
    }

    ds["flags"].encoding = {
        "dtype": "uint8",
        "zlib": True,
    }

    ds["x"].encoding = {
        "dtype": "float64",
    }

    ds["y"].encoding = {
        "dtype": "float64",
    }

    ds["time"].encoding = {
        "units": "days since 2000-01-01",
        "calendar": "standard",
    }

    return ds

@pytest.fixture
def spatial_dataset():
    """Create a dataset containing continuous and discrete spatial variables."""

    return xr.Dataset(
        data_vars={
            "reflectance": (
                ("y", "x"),
                np.array(
                    [
                        [1.0, 3.0, 5.0, 7.0],
                        [5.0, 7.0, 9.0, 11.0],
                        [9.0, 11.0, 13.0, 15.0],
                        [13.0, 15.0, 17.0, 19.0],
                    ],
                    dtype=np.float32,
                ),
            ),
            "temperature": (
                ("y", "x"),
                np.array(
                    [
                        [10.0, np.nan, 30.0, 40.0],
                        [30.0, 50.0, 50.0, 60.0],
                        [70.0, 90.0, 100.0, 110.0],
                        [90.0, 110.0, 120.0, 130.0],
                    ],
                    dtype=np.float32,
                ),
            ),
            "flags": (
                ("y", "x"),
                np.array(
                    [
                        [1, 2, 3, 4],
                        [5, 6, 7, 8],
                        [9, 10, 11, 12],
                        [13, 14, 15, 16],
                    ],
                    dtype=np.uint8,
                ),
            ),
            "mask": (
                ("y", "x"),
                np.array(
                    [
                        [True, False, True, False],
                        [False, True, False, True],
                        [True, True, False, False],
                        [False, False, True, True],
                    ],
                    dtype=bool,
                ),
            ),
        },
        coords={
            "y": [3.5, 2.5, 1.5, 0.5],
            "x": [0.5, 1.5, 2.5, 3.5],
        },
        attrs={
            "satellite": "S2A",
            "description": "Synthetic test dataset",
        },
    )


@pytest.fixture
def odd_spatial_dataset():
    """Create a spatial dataset with dimensions that are not divisible by the reduction factor."""

    return xr.Dataset(
        data_vars={
            "reflectance": (
                ("y", "x"),
                np.arange(15, dtype=np.float32).reshape(3, 5),
            ),
            "flags": (
                ("y", "x"),
                np.arange(15, dtype=np.uint8).reshape(3, 5),
            ),
        },
        coords={
            "y": [2.0, 1.0, 0.0],
            "x": [0.0, 1.0, 2.0, 3.0, 4.0],
        },
        attrs={"description": "Odd-sized dataset"},
    )
    
# ============================================================
# constants
# ============================================================

class TestConstants:
    
    def test_format_constants(self): 
        """Vérifie les constantes relatives aux formats de sortie."""

        assert NETCDF_ENGINE == "netcdf4"

        assert _CONFLICTING_ATTRS == (
            "grid_mapping",
            "scale_factor",
            "_FillValue",
            "add_offset",
        )

        assert _NETCDF_EXTENSIONS == (".nc", ".nc4", ".cdf")
        assert _ZARR_EXTENSIONS == (".zarr",)

        assert _FORMAT_EXTENSION == {
            "netcdf": ".nc",
            "zarr": ".zarr",
        }

    def test_pyramid_constants(self): 
        """Vérifie les paramètres de génération des pyramides."""

        assert PYRAMID_FACTOR == 2
        assert PYRAMID_MIN_SIZE == 512

        assert LEVEL_GROUP == "{:d}"
        assert ANCILLARY_GROUP == "ancillary"

    def test_zarr_layout_constants(self): 
        """Vérifie les paramètres de shards et de chunks Zarr."""

        assert _DEFAULT_SHARDS == {
            "y": 1024,
            "x": 1024,
            "wl": -1,
        }

        assert _DEFAULT_CHUNKS_SHARDED == {
            "y": 256,
            "x": 256,
            "wl": -1,
        }

        assert _DEFAULT_CHUNKS == {
            "y": 512,
            "x": 512,
            "wl": -1,
        }

# ============================================================
# Test Base encoding
# ============================================================

class TestBaseEncoding:
    
    def test_content(self): 
        """Check the encoding of main variables"""

        encoding = _base_encoding()

        assert encoding["aot550"] == {
            "dtype": "int16",
            "scale_factor": 0.001,
            "_FillValue": -9999,
            "grid_mapping": "spatial_ref",
        }

        assert encoding["BRDFg"] == {
            "dtype": "int16",
            "scale_factor": 0.00001,
            "add_offset": 0.3,
            "_FillValue": -32768,
            "grid_mapping": "spatial_ref",
        }

        assert encoding["dem"] == {
            "dtype": "int16",
            "scale_factor": 1,
            "add_offset": 500,
            "_FillValue": -32768,
            "grid_mapping": "spatial_ref",
        }

        assert encoding["vza"] == {
            "dtype": "int16",
            "scale_factor": 0.01,
            "_FillValue": -32768,
            "grid_mapping": "spatial_ref",
        }

        assert encoding["raa"] == {
            "dtype": "int16",
            "scale_factor": 0.01,
            "add_offset": 180,
            "_FillValue": -32768,
            "grid_mapping": "spatial_ref",
        }

        assert encoding["sza"] == {
            "dtype": "int16",
            "scale_factor": 0.01,
            "add_offset": 30,
            "_FillValue": -32768,
            "grid_mapping": "spatial_ref",
        }

        assert encoding["flags"] == {
            "grid_mapping": "spatial_ref",
        }

        assert encoding["mask"] == {
            "grid_mapping": "spatial_ref",
        }

    def test_returns_independent_dictionaries(self): 
        """Check that too calls do not share their dictionaries"""

        encoding1 = _base_encoding()
        encoding2 = _base_encoding()

        encoding1["aot550"]["scale_factor"] = 123

        assert encoding2["aot550"]["scale_factor"] == 0.001

    def test_variables(self): 
        """Check that the exact expected variables are there."""

        encoding = _base_encoding()

        assert set(encoding) == {
            "aot550",
            "BRDFg",
            "dem",
            "vza",
            "raa",
            "sza",
            "wind",
            "flags",
            "mask",
        }

# ============================================================
# Test Rrs encoding
# ============================================================

class TestRrsEncoding:
    
    def test_content(self): 
        """Vérifie l'encodage des réflectances Rrs."""

        encoding = _rrs_encoding()

        assert encoding == {
            "dtype": "int16",
            "scale_factor": 0.00001,
            "add_offset": 0.3,
            "_FillValue": -32768,
            "grid_mapping": "spatial_ref",
        }

    def test_returns_independent_dictionaries(self): 
        """Vérifie que chaque appel renvoie un nouveau dictionnaire."""

        encoding1 = _rrs_encoding()
        encoding2 = _rrs_encoding()

        encoding1["dtype"] = "float32"

        assert encoding2["dtype"] == "int16"

# ============================================================
# Test Coord encoding
# ============================================================

class TestCoordEncoding:
   
    def test_content(self): 
        """Vérifie l'encodage des coordonnées temporelles et spatiales."""

        encoding = _coord_encoding()

        assert encoding == {
            "time": {
                "units": "seconds since 1970-01-01T00:00:00",
                "calendar": "proleptic_gregorian",
                "dtype": "int64",
            },
            "spatial_ref": {
                "dtype": "int32",
            },
        }

    def test_returns_independent_dictionaries(self):
        """Vérifie l'indépendance des dictionnaires de coordonnées."""

        encoding1 = _coord_encoding()
        encoding2 = _coord_encoding()

        encoding1["time"]["dtype"] = "float64"

        assert encoding2["time"]["dtype"] == "int64"

# ============================================================
# Test _package_version
# ============================================================

class TestPackageVersion:
   
    def test_returns_installed_version(self): 
        """Returns the package version when it is installed"""

        with patch(
            "importlib.metadata.version",
            return_value="1.2.3",
        ) as mock_version:

            result = _package_version("mon-paquet")

        assert result == "1.2.3"
        mock_version.assert_called_once_with("mon-paquet")

    def test_returns_unknown_if_not_installed(self): 
        """Returns 'unknown' if the package is not installed."""

        from importlib.metadata import PackageNotFoundError

        with patch(
            "importlib.metadata.version",
            side_effect=PackageNotFoundError,
        ):

            result = _package_version("paquet-inexistant")

        assert result == "unknown"

# ============================================================
# Test _reset_encodingg
# ============================================================

class TestResetEncoding:
   
    def test_clears_all_variable_encodings(self,sample_dataset):
        """Check all encodings are removed."""

        _reset_encoding(sample_dataset)

        for var in sample_dataset.variables:
            assert sample_dataset[var].encoding == {}
        

    def test_preserves_data(self,sample_dataset): 
        """Check the data stays unchanged."""

        original = sample_dataset.copy(deep=True)

        _reset_encoding(sample_dataset)

        xr.testing.assert_equal(sample_dataset, original)

    def test_reset_encoding_returns_none(self,sample_dataset): 
        """Vérifie que la fonction modifie le Dataset en place."""

        result = _reset_encoding(sample_dataset)

        assert result is None

    def test_handles_dataset_without_encoding(self): 
        """Check the behaviour with a Dataset without any explicit encoding"""

        ds = xr.Dataset(
            {
                "temperature": (
                    ("x",),
                    np.array([10.0, 20.0]),
                )
            }
        )

        _reset_encoding(ds)

        for var in ds.variables:
            assert ds[var].encoding == {}

# ============================================================
# Test _load_coords
# ============================================================

class TestLoadCoords:
   
    def test_without_lazy_coordinates(self): 
        """Returns the original Dataset if no scalar coordinate is lazy."""

        ds = xr.Dataset(
            {"temperature": (("y", "x"), [[10.0, 20.0], [30.0, 40.0]])},
            coords={
                "x": [0, 1],
                "y": [0, 1],
                "time": np.datetime64("2024-06-15"),
            },
        )

        result = _load_coords(ds)

        assert result is ds

    def test_computes_lazy_scalar_coordinates(self):
        """Test that lazy scalar coordinates are computed and loaded into memory."""

        da = pytest.importorskip("dask.array")

        time_value = np.datetime64("2024-06-15", "ns")

        ds = xr.Dataset(
            {"temperature": (("y", "x"), np.ones((2, 2)))},
            coords={
                "x": [0, 1],
                "y": [0, 1],
                "time": xr.Variable(
                    (),
                    da.from_array(time_value, chunks=()),
                ),
            },
        )

        assert ds["time"].chunks is not None

        result = mod._load_coords(ds)

        assert result is not ds
        assert result["time"].chunks is None
        assert result["time"].shape == ()
        assert result["time"].dtype == np.dtype("datetime64[ns]")
        assert result["time"].values == time_value

        xr.testing.assert_equal(
            result["temperature"],
            ds["temperature"],
        )


    def test_does_not_compute_dimension_coordinates(self):
        """Test that dimension coordinates are not computed."""

        da = pytest.importorskip("dask.array")

        ds = xr.Dataset(
            {"temperature": (("y", "x"), np.ones((2, 2)))},
            coords={
                "x": [0, 1],
                "y": [0, 1],
            },
        )

        mocked_coord = MagicMock()
        mocked_coord.chunks = ((1, 1),)

        with patch.object(
            xr.Dataset,
            "coords",
            new_callable=lambda: property(
                lambda self: {"x": mocked_coord}
            ),
        ):
            result = mod._load_coords(ds)

        assert result is ds
        mocked_coord.compute.assert_not_called()


# ============================================================
# Test _format_from_name
# ============================================================

@pytest.mark.parametrize( ("path", "expected"), [ ("product.nc", "netcdf"), ("product.NC", "netcdf"), ("product.nc4", "netcdf"), ("product.cdf", "netcdf"), ("product.zarr", "zarr"), ("product.ZARR", "zarr"), ("/tmp/product.nc", "netcdf"), ("/tmp/product.zarr/", "zarr"), ("s3://bucket/product.nc", "netcdf"), ("s3://bucket/product.zarr/", "zarr"), ("product.tif", None), ("product", None), ("", None), ], ) 
def test_format_from_name(path, expected): 
    """Détecte correctement le format selon l'extension."""
    assert _format_from_name(path) == expected

def test_format_from_name_accepts_pathlike(): 
    """Accepte un objet pathlib.Path."""

    from pathlib import Path

    assert _format_from_name(Path("product.nc")) == "netcdf"

# ============================================================
# Test _is_remote
# ============================================================
   
@pytest.mark.parametrize( ("path", "expected"), [ ("product.nc", False), ("/tmp/product.zarr", False), ("./product.zarr", False), ("s3://bucket/product.zarr", True), ("gs://bucket/product.zarr", True), ("az://container/product.zarr", True), ("https://example.org/product.zarr", True), ("file:///tmp/product.zarr", False), ("", False), ], ) 
def test_is_remote(path, expected): 
    """Distingue les URL distantes des chemins locaux."""
    assert _is_remote(path) is expected

def test_is_remote_accepts_pathlike(): 
    """Accepte un objet pathlib.Path."""

    from pathlib import Path

    assert _is_remote(Path("/tmp/product.zarr")) is False

# ============================================================
# Test _zarr_major
# ============================================================

class TestZarrMajor:
    
    @pytest.mark.parametrize( ("version", "expected"), [ ("2.18.7", 2), ("3.0.1", 3), ("10.2.0", 10), ], ) 
    def test_returns_major_version(self,version, expected): 
        """Extract the major number of the zarr's version"""
        """Extrait le numéro majeur de la version de Zarr."""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.__version__ = version

        with patch.dict(sys.modules, {"zarr": fake_zarr}):
            assert _zarr_major() == expected
        
    def test_returns_none_if_zarr_is_missing(self): 
        """Return None if Zarr cannot be imported"""

        original_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "zarr":
                raise ImportError("zarr not installed")
            return original_import(name, *args, **kwargs)

        with (
            patch.dict(sys.modules, {"zarr": None}),
            patch("builtins.__import__", side_effect=fake_import),
        ):
            assert _zarr_major() is None

    def test_if_version_is_invalid(self): 
        """Returns None if the version cannot be interpreted"""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.__version__ = "invalid"

        with patch.dict(sys.modules, {"zarr": fake_zarr}):
            assert _zarr_major() is None

# ============================================================
# Test _xarray_supports_shards
# ============================================================
 
class TestXarraySupportsShards:
        
    @pytest.mark.parametrize( ("source", "expected"), [ ("def encode():\n shards = None", True), ("def encode():\n chunks = None", False), ], ) 
    def test_nominal(self,source, expected): 
        """Detect the presence of shards functionnality in source code"""

        with patch("inspect.getsource", return_value=source):
            assert _xarray_supports_shards() is expected

    def test_returns_false_on_error(self): 
        """Return False if code inspection fails."""

        with patch(
            "inspect.getsource",
            side_effect=OSError("source unavailable"),
        ):
            assert _xarray_supports_shards() is False


# ============================================================
# Test _compression_encoding
# ============================================================
 
class TestCompressionEncoding:
        
    @pytest.mark.parametrize( ("compression", "clevel", "zarr_major"), [ (None, 5, 2), ("zstd", 5, None), ], ) 
    def test_disabled(self,compression, clevel, zarr_major): 
        """Does not configure compression if it is deactivated or unavailable."""

        assert _compression_encoding(
            compression, clevel, zarr_major
        ) == {}

    def test_zarr_v3(self): 
        """Configure BloscCodec with Zarr 3 API."""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.__path__ = []

        fake_codecs = types.ModuleType("zarr.codecs")
        fake_codec = MagicMock(name="BloscCodec")
        fake_codecs.BloscCodec = fake_codec

        with patch.dict(
            sys.modules,
            {
                "zarr": fake_zarr,
                "zarr.codecs": fake_codecs,
            },
        ):
            result = _compression_encoding("zstd", 5, 3)

        fake_codec.assert_called_once_with(
            cname="zstd",
            clevel=5,
            shuffle="shuffle",
        )

        assert result == {
            "compressors": (fake_codec.return_value,),
        }

    def test_zarr_v2(self): 
        """Configure Blosc with  Zarr 2 API."""

        fake_numcodecs = types.ModuleType("numcodecs")
        fake_numcodecs.Blosc = MagicMock(name="Blosc")
        fake_numcodecs.Blosc.SHUFFLE = 1

        with patch.dict(sys.modules, {"numcodecs": fake_numcodecs}):
            result = _compression_encoding("zstd", 7, 2)

        fake_numcodecs.Blosc.assert_called_once_with(
            cname="zstd",
            clevel=7,
            shuffle=1,
        )

        assert result == {
            "compressor": fake_numcodecs.Blosc.return_value,
        }

    def test_returns_empty_dict_on_error(self,caplog): 
        """Log a warning and returns {} if compression fails"""

        with (
            patch.dict(sys.modules, {"numcodecs": None}),
            caplog.at_level(logging.WARNING),
        ):
            result = _compression_encoding("invalid", 5, 2)

        assert result == {}
        assert "compression 'invalid' not available" in caplog.text

# ============================================================
# Test _var_chunking
# ============================================================

class TestVarChunking:
    
    @pytest.mark.parametrize( ("dims", "dask_chunks", "chunks", "expected"), [ ( ("y", "x"), ((512, 512), (512, 512)), {"y": 256, "x": 256}, ((256, 256), (512, 512)), ), ( ("y", "x", "wl"), ((512, 512), (512, 512), (10,)), {"y": 256, "x": 256, "wl": -1}, ((256, 256, 10), (512, 512, 10)), ), ( ("y", "x"), ((100, 100), (120, 120)), {"y": 256, "x": 256}, ((100, 120), (100, 120)), ), ( ("y", "x"), ((500, 500), (300, 300)), {"y": 256, "x": 100}, ((500, 100), (500, 300)), ), ( ("y", "x"), ((512, 512), (512, 512)), {}, ((512, 512), (512, 512)), ), ( ("y", "x"), ((512, 512), (512, 512)), {"y": None, "x": -1}, ((512, 512), (512, 512)), ), ], ) 
    def test_nominal(self,dims, dask_chunks, chunks, expected): 
        """Checks the calculation of stored chunks"""

        assert _var_chunking(dims, dask_chunks, chunks) == expected

    def test_clips_chunk_to_block(self): 
        """Limit the chunk to the size of the block."""

        result = _var_chunking(
            dims=("y",),
            dask_chunks=((128, 128),),
            chunks={"y": 512},
        )

        assert result == ((128,), (128,))

# ============================================================
# Test _clear_store
# ============================================================

class TestClearStore:
    
    def test_removes_existing_local_store(self,tmp_path):
        """Deletes an existing local store."""

        store = tmp_path / "product.zarr"
        store.mkdir()
        (store / "metadata.txt").write_text("test")

        _clear_store(str(store))

        assert not store.exists()

    def test_does_nothing_if_local_store_missing(self,tmp_path): 
        """Does nothing if local store does not exist."""

        store = tmp_path / "missing.zarr"

        with patch.object(shutil, "rmtree") as mock_rmtree:
            _clear_store(str(store))

        mock_rmtree.assert_not_called()

    def test_removes_existing_remote_store(self): 
        """Deletes a distant store via fsspec."""

        fake_fsspec = types.ModuleType("fsspec")
        fake_fsspec.core = MagicMock()

        fake_fs = MagicMock()
        fake_fs.exists.return_value = True

        fake_fsspec.core.url_to_fs.return_value = (
            fake_fs,
            "bucket/product.zarr",
        )

        with patch.dict(sys.modules, {"fsspec": fake_fsspec}):
            _clear_store(
                "s3://bucket/product.zarr",
                storage_options={"anon": True},
            )

        fake_fsspec.core.url_to_fs.assert_called_once_with(
            "s3://bucket/product.zarr",
            anon=True,
        )

        fake_fs.rm.assert_called_once_with(
            "bucket/product.zarr",
            recursive=True,
        )

    def test_does_not_remove_missing_remote_store(self):
        """Removes nothing if the distant store does not exist."""

        fake_fsspec = types.ModuleType("fsspec")
        fake_fsspec.core = MagicMock()

        fake_fs = MagicMock()
        fake_fs.exists.return_value = False

        fake_fsspec.core.url_to_fs.return_value = (
            fake_fs,
            "bucket/missing.zarr",
        )

        with patch.dict(sys.modules, {"fsspec": fake_fsspec}):
            _clear_store("s3://bucket/missing.zarr")

        fake_fs.rm.assert_not_called()

# ============================================================
# Test _consolidate
# ============================================================

class TestConsolidate:
    
    def test_local_store(self): 
        """Consolidates the root and groups of a local store."""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.consolidate_metadata = MagicMock()

        with patch.dict(sys.modules, {"zarr": fake_zarr}):
            _consolidate(
                "/tmp/product.zarr",
                groups=["0", "1", "ancillary"],
            )

        calls = fake_zarr.consolidate_metadata.call_args_list

        assert len(calls) == 4
        assert calls[0].args == ("/tmp/product.zarr",)
        assert calls[1].kwargs == {"path": "0"}
        assert calls[2].kwargs == {"path": "1"}
        assert calls[3].kwargs == {"path": "ancillary"}

    def test_consolidate_remote_store_zarr_v2(self): 
        """Use fsspec.get_mapper for a Zarr 2 distant store."""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.consolidate_metadata = MagicMock()

        fake_fsspec = types.ModuleType("fsspec")
        fake_fsspec.get_mapper = MagicMock(return_value="fake-mapper")

        with patch.dict(
            sys.modules,
            {
                "zarr": fake_zarr,
                "fsspec": fake_fsspec,
            },
        ):
            _consolidate(
                "s3://bucket/product.zarr",
                storage_options={"anon": True},
                zarr_major=2,
                groups=["0"],
            )

        fake_fsspec.get_mapper.assert_called_once_with(
            "s3://bucket/product.zarr",
            anon=True,
        )

        assert fake_zarr.consolidate_metadata.call_count == 2
        assert fake_zarr.consolidate_metadata.call_args_list[0].args == (
            "fake-mapper",
        )
        assert fake_zarr.consolidate_metadata.call_args_list[1].kwargs == {
            "path": "0",
        }

    def test_remote_store_zarr_v3(self):
        """Use FsspecStore.from_url for a Zarr 3 distant store."""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.consolidate_metadata = MagicMock()

        fake_zarr.storage = MagicMock()
        fake_zarr.storage.FsspecStore.from_url.return_value = "fake-store"

        with patch.dict(sys.modules, {"zarr": fake_zarr}):
            _consolidate(
                "s3://bucket/product.zarr",
                storage_options={"anon": True},
                zarr_major=3,
                groups=["ancillary"],
            )

        fake_zarr.storage.FsspecStore.from_url.assert_called_once_with(
            "s3://bucket/product.zarr",
            storage_options={"anon": True},
        )

        assert fake_zarr.consolidate_metadata.call_count == 2
        assert fake_zarr.consolidate_metadata.call_args_list[0].args == (
            "fake-store",
        )
        assert fake_zarr.consolidate_metadata.call_args_list[1].kwargs == {
            "path": "ancillary",
        }

    def test_remote_store_without_storage_options(self): 
        """Use empty options if no options is given."""

        fake_zarr = types.ModuleType("zarr")
        fake_zarr.consolidate_metadata = MagicMock()

        fake_fsspec = types.ModuleType("fsspec")
        fake_fsspec.get_mapper = MagicMock(return_value="fake-mapper")

        with patch.dict(
            sys.modules,
            {
                "zarr": fake_zarr,
                "fsspec": fake_fsspec,
            },
        ):
            _consolidate(
                "s3://bucket/product.zarr",
                zarr_major=2,
            )

        fake_fsspec.get_mapper.assert_called_once_with(
            "s3://bucket/product.zarr",
        )

# ============================================================
# Test _pyramid_levels
# ============================================================

class TestPyramidLevels:
 
    @pytest.mark.parametrize(
        ("ny", "nx", "min_size", "max_levels", "expected"),
        [
            (5490, 5490, 512, None, 4),
            (400, 300, 512, None, 0),
            (1024, 1024, 512, None, 1),
            (1025, 1025, 512, None, 2),
            (2048, 1024, 512, None, 2),
            (4096, 2048, 512, None, 3),
            (5490, 5490, 512, 2, 2),
            (5490, 5490, 512, 1, 1),
            (5490, 5490, 512, 10, 4),
            (1, 1000, 512, None, 0),
            (1000, 1, 512, None, 0),
            (0, 1000, 512, None, 0),
            (1000, 0, 512, None, 0),
            (512, 512, 512, None, 0),
            (513, 512, 512, None, 1),
        ],
    )
    def test_nominal(self,ny, nx, min_size, max_levels, expected):
        """Test the number of pyramid overview levels."""

        assert _pyramid_levels(
            ny,
            nx,
            min_size=min_size,
            max_levels=max_levels,
        ) == expected


    def test_uses_default_min_size(self):
        """Test the default minimum pyramid size."""

        assert _pyramid_levels(400, 300) == 0
        assert _pyramid_levels(1024, 1024) == 1


    def test_stops_at_max_levels(self):
        """Test that the maximum number of levels is respected."""

        assert _pyramid_levels(
            10000,
            10000,
            min_size=512,
            max_levels=2,
        ) == 2

# ============================================================
# Test _coarsen
# ============================================================

class TestCoarsen:

    def test_reduces_spatial_dimensions(self,spatial_dataset):
        """Test that spatial dimensions are reduced by the specified factor."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset, factor=2)

        assert result.sizes["y"] == 2
        assert result.sizes["x"] == 2


    def test_computes_mean_for_continuous_variables(self,spatial_dataset):
        """Test that continuous variables are averaged over each block."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset, factor=2)

        expected = np.array(
            [
                [4.0, 8.0],
                [12.0, 16.0],
            ],
            dtype=np.float32,
        )

        np.testing.assert_allclose(
            result["reflectance"].values,
            expected,
        )


    def test_ignores_nan_when_computing_mean(self,spatial_dataset):
        """Test that NaN values are ignored when averaging continuous variables."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset, factor=2)

        expected = np.array(
            [
                [30.0, 45.0],
                [90.0, 115.0],
            ],
            dtype=np.float32,
        )

        np.testing.assert_allclose(
            result["temperature"].values,
            expected,
        )


    def test_preserves_first_pixel_for_integer_variables(self,spatial_dataset):
        """Test that integer variables use the first pixel of each block."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset, factor=2)

        expected = np.array(
            [
                [1, 3],
                [9, 11],
            ],
            dtype=np.uint8,
        )

        np.testing.assert_array_equal(
            result["flags"].values,
            expected,
        )


    def test_preserves_first_pixel_for_boolean_variables(self,spatial_dataset):
        """Test that boolean variables use the first pixel of each block."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset, factor=2)

        expected = np.array(
            [
                [True, True],
                [True, False],
            ],
            dtype=bool,
        )

        np.testing.assert_array_equal(
            result["mask"].values,
            expected,
        )


    def test_preserves_dataset_attributes(self,spatial_dataset):
        """Test that dataset attributes are preserved."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset, factor=2)

        assert result.attrs == spatial_dataset.attrs


    def test_preserves_variable_attributes(self):
        """Test that continuous variable attributes are preserved."""

        ds = xr.Dataset(
            {
                "reflectance": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                    {"units": "sr-1", "description": "Remote sensing reflectance"},
                )
            },
            coords={
                "y": [3.5, 2.5, 1.5, 0.5],
                "x": [0.5, 1.5, 2.5, 3.5],
            },
        )

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(ds, factor=2)

        assert result["reflectance"].attrs == {
            "units": "sr-1",
            "description": "Remote sensing reflectance",
        }


    def test_handles_incomplete_edge_blocks(self,odd_spatial_dataset):
        """Test that incomplete blocks at the right and bottom edges are retained."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(odd_spatial_dataset, factor=2)

        assert result.sizes["y"] == 2
        assert result.sizes["x"] == 3

        np.testing.assert_allclose(
            result["reflectance"].values,
            np.array(
                [
                    [3.0, 5.0, 6.5],
                    [10.5, 12.5, 14.0],
                ],
                dtype=np.float32,
            ),
        )

        np.testing.assert_array_equal(
            result["flags"].values,
            np.array(
                [
                    [0, 2, 4],
                    [10, 12, 14],
                ],
                dtype=np.uint8,
            ),
        )


    def test_updates_geotransform(self,odd_spatial_dataset):
        """Test that the geotransform update function is called."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ) as mock_update:
            result = _coarsen(odd_spatial_dataset, factor=2)

        mock_update.assert_called_once()

        assert result.sizes["y"] == 2
        assert result.sizes["x"] == 3


    def test_uses_default_factor(self,spatial_dataset):
        """Test that the default reduction factor is applied."""

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(spatial_dataset)

        assert result.sizes["y"] == 2
        assert result.sizes["x"] == 2


    def test_preserves_non_spatial_variables(self):
        """Test the behavior of variables that do not have both spatial dimensions."""

        ds = xr.Dataset(
            data_vars={
                "reflectance": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                ),
                "quality": (
                    ("wl",),
                    np.array([1, 2, 3], dtype=np.int16),
                ),
            },
            coords={
                "y": [3.5, 2.5, 1.5, 0.5],
                "x": [0.5, 1.5, 2.5, 3.5],
                "wl": [490, 560, 842],
            },
        )

        with patch.object(
            mod,
            "_update_geotransform",
            side_effect=lambda ds: ds,
        ):
            result = _coarsen(ds, factor=2)

        assert "quality" in result
        assert result["quality"].sizes == {"wl": 3}

        np.testing.assert_array_equal(
            result["quality"].values,
            [1, 2, 3],
        )


# ============================================================
# Test _update_geotransform
# ============================================================

class TestUpdateGeotransform:

    def test_updates_transform(self):
        """Test that spatial dimensions and geotransform are updated."""

        ds = xr.Dataset(
            {
                "reflectance": (
                    ("row", "col"),
                    np.ones((2, 3), dtype=np.float32),
                )
            },
            coords={
                "row": [1.5, 0.5],
                "col": [0.5, 1.5, 2.5],
            },
        )

        expected_transform = MagicMock(name="transform")

        rio_mock = MagicMock()
        rio_mock.transform.return_value = expected_transform
        rio_mock.set_spatial_dims.return_value = ds
        rio_mock.write_transform.return_value = ds

        with patch.object(
            xr.Dataset,
            "rio",
            new_callable=lambda: property(lambda self: rio_mock),
            create=True,
        ):
            result = _update_geotransform(
                ds,
                dims=("row", "col"),
            )

        rio_mock.set_spatial_dims.assert_called_once_with(
            x_dim="col",
            y_dim="row",
        )
        rio_mock.transform.assert_called_once_with(recalc=True)
        rio_mock.write_transform.assert_called_once_with(expected_transform)

        assert result is ds


    def test_uses_default_dimensions(self):
        """Test that the default spatial dimension names are used."""

        ds = xr.Dataset(
            {
                "reflectance": (
                    ("y", "x"),
                    np.ones((2, 2), dtype=np.float32),
                )
            },
            coords={
                "y": [1.5, 0.5],
                "x": [0.5, 1.5],
            },
        )

        rio_mock = MagicMock()
        rio_mock.set_spatial_dims.return_value = ds
        rio_mock.transform.return_value = "fake-transform"
        rio_mock.write_transform.return_value = ds

        with patch.object(
            xr.Dataset,
            "rio",
            new_callable=lambda: property(lambda self: rio_mock),
            create=True,
        ):
            result = _update_geotransform(ds)

        rio_mock.set_spatial_dims.assert_called_once_with(
            x_dim="x",
            y_dim="y",
        )

        assert result is ds


    def test_returns_dataset_if_spatial_dims_fail(self,caplog):
        """Test that the original dataset is returned if spatial dimension setup fails."""

        ds = xr.Dataset(
            {
                "reflectance": (
                    ("y", "x"),
                    np.ones((2, 2), dtype=np.float32),
                )
            }
        )

        rio_mock = MagicMock()
        rio_mock.set_spatial_dims.side_effect = ValueError("invalid spatial dimensions")

        with (
            patch.object(
                xr.Dataset,
                "rio",
                new_callable=lambda: property(lambda self: rio_mock),
                create=True,
            ),
            caplog.at_level(logging.WARNING),
        ):
            result = _update_geotransform(ds)

        assert result is ds
        assert "could not update the geotransform" in caplog.text
        assert "invalid spatial dimensions" in caplog.text


    def test_returns_dataset_if_transform_fails(self,caplog):
        """Test that the dataset is returned if calculating the transform fails."""

        ds = xr.Dataset(
            {
                "reflectance": (
                    ("y", "x"),
                    np.ones((2, 2), dtype=np.float32),
                )
            }
        )

        rio_mock = MagicMock()
        rio_mock.set_spatial_dims.return_value = ds
        rio_mock.transform.side_effect = RuntimeError("transform unavailable")

        with (
            patch.object(
                xr.Dataset,
                "rio",
                new_callable=lambda: property(lambda self: rio_mock),
                create=True,
            ),
            caplog.at_level(logging.WARNING),
        ):
            result = _update_geotransform(ds)

        assert result is ds
        assert "could not update the geotransform" in caplog.text
        assert "transform unavailable" in caplog.text


    def test_returns_dataset_if_write_transform_fails(self,caplog):
        """Test that the dataset is returned if writing the transform fails."""

        ds = xr.Dataset(
            {
                "reflectance": (
                    ("y", "x"),
                    np.ones((2, 2), dtype=np.float32),
                )
            }
        )

        rio_mock = MagicMock()
        rio_mock.set_spatial_dims.return_value = ds
        rio_mock.transform.return_value = "fake-transform"
        rio_mock.write_transform.side_effect = RuntimeError("write failed")

        with (
            patch.object(
                xr.Dataset,
                "rio",
                new_callable=lambda: property(lambda self: rio_mock),
                create=True,
            ),
            patch.object(
                logging,
                "warning",
            ) as mock_warning,
        ):
            result = _update_geotransform(ds)

        assert result is ds
        mock_warning.assert_called_once()
        assert "could not update the geotransform" in mock_warning.call_args.args[0]