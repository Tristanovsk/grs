
import numpy as np 
import pytest 
import xarray as xr

import grs.output as mod
from unittest.mock import MagicMock
from grs.output import L2aProduct

# ============================================================
# Fixtures and helpers
# ============================================================

class FakeRioAccessor:
    """Minimal replacement for the rioxarray accessor used by unit tests."""

    def __init__(self, obj):
        self._obj = obj

    def set_spatial_dims(self, x_dim="x", y_dim="y", inplace=False):
        self._obj.attrs["_x_dim"] = x_dim
        self._obj.attrs["_y_dim"] = y_dim
        return self._obj

    def write_coordinate_system(self, inplace=False):
        return self._obj

    def write_crs(self, input_crs=None, grid_mapping_name="spatial_ref", inplace=False, **kwargs):
        if input_crs is not None:
            self._obj.attrs["crs"] = str(input_crs)
        return self._obj

    def transform(self, recalc=False):
        return (0, 1, 0, 0, 0, -1)

    def write_transform(self, transform=None, grid_mapping_name="spatial_ref", inplace=False):
        if transform is not None:
            self._obj.attrs["GeoTransform"] = str(transform)
        return self._obj


@pytest.fixture(autouse=True)
def fake_rio_accessor(monkeypatch):
    """Register a minimal rio accessor without requiring rioxarray."""

    def get_rio(obj):
        return FakeRioAccessor(obj)

    monkeypatch.setattr(
        xr.Dataset,
        "rio",
        property(get_rio),
        raising=False,
    )
    monkeypatch.setattr(
        xr.DataArray,
        "rio",
        property(get_rio),
        raising=False,
    )


@pytest.fixture(autouse=True)
def disable_geotransform_update(monkeypatch):
    """Disable geotransform updates that depend on rioxarray."""

    monkeypatch.setattr(
        "grs.output._update_geotransform",
        lambda ds, dims=("y", "x"): ds,
    )


@pytest.fixture
def sample_products():
    """Create synthetic products for testing L2aProduct."""

    x = np.array([100.0, 110.0, 120.0, 130.0])
    y = np.array([200.0, 190.0, 180.0])
    wl = np.array([490, 560, 665])

    native_raster = xr.Dataset(
        {
            "vza": (
                ("wl", "y", "x"),
                np.full((3, 3, 4), 10.0, dtype=np.float32),
            ),
            "sza": (
                ("y", "x"),
                np.full((3, 4), 30.0, dtype=np.float32),
            ),
            "raa": (
                ("wl", "y", "x"),
                np.full((3, 3, 4), 20.0, dtype=np.float32),
            ),
            "flags": (
                ("y", "x"),
                np.zeros((3, 4), dtype=np.uint8),
            ),
            "wl_true": (("wl",), wl),
        },
        coords={"x": x, "y": y, "wl": wl},
        attrs={
            "satellite": "S2A",
            "acquisition_date": "2024-06-15T10:30:00.000Z",
            "test_attribute": "native_metadata",
        },
    )

    l2_prod = xr.Dataset(
        {
            "Rrs": (
                ("wl", "y", "x"),
                np.full((3, 3, 4), 0.01, dtype=np.float32),
            ),
            "aot550": (
                ("y", "x"),
                np.full((3, 4), 0.1, dtype=np.float32),
            ),
            "BRDFg": (
                ("y", "x"),
                np.full((3, 4), 0.3, dtype=np.float32),
            ),
        },
        coords={"x": x, "y": y, "wl": wl},
    )

    cams_raster = xr.Dataset(
        {
            "aot_cams": (
                ("y", "x"),
                np.full((2, 2), 0.2, dtype=np.float32),
            ),
        },
        coords={
            "x": [105.0, 125.0],
            "y": [195.0, 175.0],
        },
    )

    gas_trans = MagicMock()
    gas_trans.Tg_tot_coarse = xr.DataArray(
        np.ones((2, 2, 3), dtype=np.float32),
        dims=("x", "y", "wl"),
        coords={
            "x": [105.0, 125.0],
            "y": [195.0, 175.0],
            "wl": wl,
        },
        name="Tg_tot_coarse",
    )

    prod = MagicMock()
    prod.raster = native_raster
    prod.sunglint_threshold = 0.2
    prod.ndwi_threshold = 0.0
    prod.vis_swir_index_threshold = 0.1
    prod.hcld_threshold = 0.003
    prod.dirdata = "/tmp/grsdata"
    prod.abs_gas_file = "/tmp/abs_gas.nc"
    prod.water_vapor_transmittance_file = "/tmp/water_vapor.nc"
    prod.bcirrus = None
    prod.cirrus = None
    prod.bwv = None
    prod.wv = None

    return prod, l2_prod, cams_raster, gas_trans


@pytest.fixture
def l2a_product(sample_products):
    """Create an initialized L2aProduct."""

    prod, l2_prod, cams_raster, gas_trans = sample_products

    return L2aProduct(
        prod=prod,
        l2_prod=l2_prod,
        cams=MagicMock(raster=cams_raster),
        gas_trans=gas_trans,
    )

# ============================================================
# L2A Product initialization and construction
# ============================================================

class TestL2AProductConstruction:
    
    def test_initialization(self,sample_products):
        """Test initialization of the L2aProduct object."""

        prod, l2_prod, cams_raster, gas_trans = sample_products

        product = L2aProduct(
            prod=prod,
            l2_prod=l2_prod,
            cams=MagicMock(raster=cams_raster),
            gas_trans=gas_trans,
        )

        assert product.prod is prod
        assert product.cams.raster is cams_raster
        assert product.gas_trans is gas_trans
        assert product.dem is None
        assert product.ancillary is not None
        assert product.wl_cirrus_band == 1375
        assert product.wl_wv_band == 945
        assert product.complevel == 5

    def test_main_product_variables(self,l2a_product):
        """Test that the main product contains the expected variables."""

        expected_variables = {
            "Rrs",
            "aot550",
            "BRDFg",
            "vza",
            "sza",
            "raa",
            "flags",
            "dem",
        }

        assert expected_variables.issubset(set(l2a_product.l2_prod.data_vars))

    def test_metadata(self,l2a_product):
        """Test that native metadata and processing parameters are copied."""

        attrs = l2a_product.l2_prod.attrs

        assert attrs["satellite"] == "S2A"
        assert attrs["acquisition_date"] == "2024-06-15T10:30:00.000Z"
        assert attrs["test_attribute"] == "native_metadata"
        assert attrs["sunglint_threshold"] == "0.2"
        assert attrs["ndwi_threshold"] == "0.0"
        assert attrs["vis_swir_index_threshold"] == "0.1"
        assert attrs["hcld_threshold"] == "0.003"
        assert attrs["dirdata"] == "/tmp/grsdata"
        assert attrs["abs_gas_file"] == "/tmp/abs_gas.nc"
        assert attrs["water_vapor_transmittance_file"] == "/tmp/water_vapor.nc"
        assert "GRSdriver_version" in attrs

    def test_does_not_modify_native_attributes(self,sample_products):
        """Test that native raster attributes remain unchanged."""

        prod, l2_prod, cams_raster, gas_trans = sample_products
        original_attrs = dict(prod.raster.attrs)

        L2aProduct(
            prod=prod,
            l2_prod=l2_prod,
            cams=MagicMock(raster=cams_raster),
            gas_trans=gas_trans,
        )

        assert prod.raster.attrs == original_attrs

    def test_extracts_angles(self,l2a_product):
        """Test that geometry angles are extracted correctly."""

        ds = l2a_product.l2_prod

        np.testing.assert_allclose(ds["vza"].values, 10.0)
        np.testing.assert_allclose(ds["sza"].values, 30.0)
        np.testing.assert_allclose(ds["raa"].values, 20.0)

        assert ds["vza"].dims == ("y", "x")
        assert ds["sza"].dims == ("y", "x")
        assert ds["raa"].dims == ("y", "x")

    def test_creates_empty_dem_when_missing(self,l2a_product):
        """Test that a missing DEM is replaced with a NaN-filled array."""

        dem = l2a_product.l2_prod["dem"]

        assert dem.dtype == np.float32
        assert dem.dims == ("y", "x")
        assert bool(dem.isnull().all())

    def test_uses_provided_dem(self,sample_products):
        """Test that the provided DEM is included in the main product."""

        prod, l2_prod, cams_raster, gas_trans = sample_products

        dem = xr.DataArray(
            np.full((3, 4), 150.0, dtype=np.float32),
            dims=("y", "x"),
            coords={"x": l2_prod.x, "y": l2_prod.y},
            name="dem",
        )

        product = L2aProduct(
            prod=prod,
            l2_prod=l2_prod,
            cams=MagicMock(raster=cams_raster),
            gas_trans=gas_trans,
            dem=dem,
        )

        np.testing.assert_array_equal(
            product.l2_prod["dem"].values,
            dem.values,
        )

    def test_ancillary_variables(self,l2a_product):
        """Test that ancillary variables use coarse-resolution dimensions."""

        ancillary = l2a_product.ancillary

        assert "aot_cams" in ancillary
        assert "Tg_tot_coarse" in ancillary
        assert "xc" in ancillary.dims
        assert "yc" in ancillary.dims
        assert "x" not in ancillary.dims
        assert "y" not in ancillary.dims
        assert "wl_all" in ancillary.dims


# ============================================================
# L2A Attributes
# ============================================================

class TestL2AProductAttributes:
    
    def test_scrub_attrs_removes_conflicting_attributes(self):
        """Test that conflicting variable attributes are removed."""

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((2, 2), dtype=np.float32),
                    {
                        "grid_mapping": "spatial_ref",
                        "scale_factor": 0.00001,
                        "_FillValue": -32768,
                        "add_offset": 0.3,
                        "units": "sr-1",
                    },
                ),
            }
        )

        L2aProduct._scrub_attrs(ds)

        for attr in (
            "grid_mapping",
            "scale_factor",
            "_FillValue",
            "add_offset",
        ):
            assert attr not in ds["Rrs"].attrs

        assert ds["Rrs"].attrs["units"] == "sr-1"


    def test_bool_attrs_to_int_converts_dataset_and_variable_attributes(self):
        """Test that boolean attributes are converted to integers."""

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((2, 2)),
                    {"valid": True, "invalid": False},
                ),
            },
            attrs={
                "flag_a": True,
                "flag_b": False,
                "description": "test",
            },
        )

        L2aProduct._bool_attrs_to_int(ds)

        assert ds.attrs["flag_a"] == 1
        assert ds.attrs["flag_b"] == 0
        assert type(ds.attrs["flag_a"]) is int
        assert ds["Rrs"].attrs["valid"] == 1
        assert ds["Rrs"].attrs["invalid"] == 0
        assert type(ds["Rrs"].attrs["valid"]) is int
        assert ds.attrs["description"] == "test"


# ============================================================
# Test _prepare_export
# ============================================================

class TestPrepareExport:
    
    def test_datacube(self,l2a_product):
        """Test preparation of the datacube export."""

        ds, encoding = l2a_product._prepare_export(snap_compliant=False)

        assert ds.attrs["metadata_profile"] == "datacube"
        assert "Rrs" in ds.data_vars
        assert "Rrs" in encoding
        assert encoding["Rrs"]["dtype"] == "int16"
        assert encoding["Rrs"]["scale_factor"] == 0.00001
        assert encoding["Rrs"]["add_offset"] == 0.3
        assert encoding["Rrs"]["_FillValue"] == -32768
        assert "aot550" in encoding
        assert "BRDFg" in encoding
        assert "grid_mapping" not in ds["Rrs"].attrs


    def test_snap_compliant(self,l2a_product):
        """Test preparation of the SNAP-compatible export."""

        ds, encoding = l2a_product._prepare_export(snap_compliant=True)

        assert ds.attrs["metadata_profile"] == "beam"
        assert ds.attrs["auto_grouping"] == "Rrs:Rrs_g"
        assert "Rrs" not in ds.data_vars

        for wavelength in (490, 560, 665):
            variable = f"Rrs_{wavelength}"

            assert variable in ds.data_vars
            assert variable in encoding
            assert encoding[variable]["dtype"] == "int16"
            assert encoding[variable]["scale_factor"] == 0.00001
            assert encoding[variable]["add_offset"] == 0.3

            attrs = ds[variable].attrs
            assert attrs["Unit"] == "sr-1"
            assert attrs["units"] == "sr-1"
            assert attrs["radiation_wavelength"] == wavelength
            assert attrs["radiation_wavelength_unit"] == "nm"
            assert attrs["wavelength"] == wavelength


    def test_does_not_modify_original_product(self,l2a_product):
        """Test that export preparation does not modify the original dataset."""

        original_attrs = dict(l2a_product.l2_prod.attrs)
        original_rrs_attrs = dict(l2a_product.l2_prod["Rrs"].attrs)

        l2a_product._prepare_export(snap_compliant=False)

        assert l2a_product.l2_prod.attrs == original_attrs
        assert l2a_product.l2_prod["Rrs"].attrs == original_rrs_attrs

# ============================================================
# Test _output_names
# ============================================================

class TestOutputNames:
    
    @pytest.mark.parametrize(
        "output_path, ext, expected_basename, expected_main",
        [
            ("scene", ".nc", "scene", "scene/scene.nc"),
            ("scene.nc", ".nc", "scene", "scene.nc"),
            ("scene.zarr", ".zarr", "scene", "scene.zarr"),
            ("scene.nc", ".zarr", "scene", "scene.nc"),
        ],
    )
    def test_nominal(
        self,
        tmp_path,
        output_path,
        ext,
        expected_basename,
        expected_main,
    ):
        """Test output path resolution and directory creation."""

        path = tmp_path / output_path

        directory, basename, main = L2aProduct._output_names(
            str(path),
            ext,
        )

        assert basename == expected_basename
        assert main == str(tmp_path / expected_main)
        assert (tmp_path / directory).is_dir()

    def test_without_extension(self,tmp_path):
        """Test legacy output path resolution without an extension."""

        output_path = tmp_path / "scene"

        directory, basename, main = L2aProduct._output_names(
            str(output_path),
            None,
        )

        assert directory == str(output_path)
        assert basename == "scene"
        assert main == str(output_path / "scene")
        assert output_path.is_dir()

# ============================================================
# Test Export configuration and verifications
# ============================================================

class TestExportConfiguration:
    
    def test_defaults_to_netcdf(self,l2a_product, monkeypatch):
        """Test that export defaults to NetCDF when no format is specified."""

        mock_export = MagicMock(return_value=None)
        monkeypatch.setattr(l2a_product, "export_to_netcdf", mock_export)

        result = l2a_product.export("scene")

        mock_export.assert_called_once_with("scene")
        assert result is None


    @pytest.mark.parametrize(
        "path, expected_format",
        [
            ("scene.nc", "netcdf"),
            ("scene.nc4", "netcdf"),
            ("scene.cdf", "netcdf"),
            ("scene.zarr", "zarr"),
            ("scene", "netcdf"),
            ("scene.tif", "netcdf"),
        ],
    )
    def test_detects_format_from_path(
        self,
        l2a_product,
        monkeypatch,
        path,
        expected_format,
    ):
        """Test that export selects the format from the path extension."""

        mock_netcdf = MagicMock(return_value="netcdf_result")
        mock_zarr = MagicMock(return_value="zarr_result")

        monkeypatch.setattr(l2a_product, "export_to_netcdf", mock_netcdf)
        monkeypatch.setattr(l2a_product, "export_to_zarr", mock_zarr)

        result = l2a_product.export(path)

        if expected_format == "netcdf":
            mock_netcdf.assert_called_once_with(path)
            mock_zarr.assert_not_called()
            assert result == "netcdf_result"
        else:
            mock_zarr.assert_called_once_with(path)
            mock_netcdf.assert_not_called()
            assert result == "zarr_result"


    @pytest.mark.parametrize("fmt", ["netcdf", "zarr"])
    def test_respects_explicit_format(self,l2a_product, monkeypatch, fmt):
        """Test that an explicit format determines the export method."""

        mock_netcdf = MagicMock(return_value=None)
        mock_zarr = MagicMock(return_value="store_path")

        monkeypatch.setattr(l2a_product, "export_to_netcdf", mock_netcdf)
        monkeypatch.setattr(l2a_product, "export_to_zarr", mock_zarr)

        result = l2a_product.export("scene", fmt=fmt, snap_compliant=True)

        if fmt == "netcdf":
            mock_netcdf.assert_called_once_with("scene", snap_compliant=True)
            mock_zarr.assert_not_called()
            assert result is None
        else:
            mock_zarr.assert_called_once_with("scene", snap_compliant=True)
            mock_netcdf.assert_not_called()
            assert result == "store_path"


    def test_passes_kwargs(self,l2a_product, monkeypatch):
        """Test that export forwards keyword arguments to the selected method."""

        mock_export = MagicMock(return_value="store_path")
        monkeypatch.setattr(l2a_product, "export_to_zarr", mock_export)

        result = l2a_product.export(
            "scene.zarr",
            chunks={"x": 128, "y": 128},
            pyramids=False,
            compression=None,
        )

        mock_export.assert_called_once_with(
            "scene.zarr",
            chunks={"x": 128, "y": 128},
            pyramids=False,
            compression=None,
        )
        assert result == "store_path"


    @pytest.mark.parametrize("fmt", ["hdf5", "nc", "ZARR", "", 123])
    def test_rejects_unknown_format(self,l2a_product, fmt):
        """Test that export rejects unsupported format names."""

        with pytest.raises(ValueError, match="unknown format"):
            l2a_product.export("scene", fmt=fmt)


    @pytest.mark.parametrize(
        "path, fmt",
        [
            ("scene.nc", "zarr"),
            ("scene.zarr", "netcdf"),
            ("scene.nc4", "zarr"),
            ("scene.cdf", "zarr"),
        ],
    )
    def test_rejects_conflicting_format_and_extension(
        self,
        l2a_product,
        path,
        fmt,
    ):
        """Test that export rejects conflicting formats and extensions."""

        with pytest.raises(ValueError, match="contradicts the extension"):
            l2a_product.export(path, fmt=fmt)


def test_export_to_netcdf_rejects_remote_paths(l2a_product):
    """Test that NetCDF export rejects remote URLs."""

    with pytest.raises(ValueError, match="NetCDF can only be written to a local path"):
        l2a_product.export_to_netcdf("s3://bucket/scene.nc")


# ============================================================
# Test Export to NetCDF
# ============================================================

class TestExportToNetCdf:
    
    def test_writes_main_and_ancillary_files(self,l2a_product, tmp_path, monkeypatch):
        calls = []

        def fake_to_netcdf(self, path, *args, **kwargs):
            calls.append({
                "path": str(path),
                "args": args,
                "kwargs": kwargs,
            })

        monkeypatch.setattr(xr.Dataset, "to_netcdf", fake_to_netcdf)

        l2a_product.export_to_netcdf(str(tmp_path / "scene.nc"))

        assert len(calls) == 2

        main_call, ancillary_call = calls

        assert main_call["path"] == str(tmp_path / "scene.nc")
        assert ancillary_call["path"] == str(tmp_path / "scene_anc.nc")

        assert main_call["kwargs"]["engine"] == mod.NETCDF_ENGINE
        assert ancillary_call["kwargs"]["engine"] == mod.NETCDF_ENGINE

        assert main_call["kwargs"]["encoding"]
        assert ancillary_call["kwargs"]["encoding"]

        assert all(
            encoding["zlib"] is True
            and encoding["complevel"] == l2a_product.complevel
            for encoding in main_call["kwargs"]["encoding"].values()
        )

        assert all(
            encoding["zlib"] is True
            and encoding["complevel"] == l2a_product.complevel
            for encoding in ancillary_call["kwargs"]["encoding"].values()
        )

    def test_adds_compression_options(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that NetCDF encodings receive the configured compression level."""

        calls = []

        def fake_to_netcdf(self, path, *args, **kwargs):
            calls.append({
                "path": str(path),
                "encoding": kwargs.get("encoding", {}),
                "engine": kwargs.get("engine"),
            })

        monkeypatch.setattr(xr.Dataset, "to_netcdf", fake_to_netcdf)
        l2a_product.complevel = 7

        l2a_product.export_to_netcdf(str(tmp_path / "scene.nc"))
        assert len(calls) == 2

        for call in calls:
            assert call["engine"] == mod.NETCDF_ENGINE

            for encoding in call["encoding"].values():
                assert encoding["zlib"] is True
                assert encoding["complevel"] == l2a_product.complevel

    def test_removes_existing_files(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that existing NetCDF outputs are removed before writing."""

        output_path = tmp_path / "scene.nc"
        ancillary_path = tmp_path / "scene_anc.nc"

        output_path.touch()
        ancillary_path.touch()

        written_files = []

        def fake_to_netcdf(self, path, *args, **kwargs):
            written_files.append(str(path))

        monkeypatch.setattr(xr.Dataset, "to_netcdf", fake_to_netcdf)

        l2a_product.export_to_netcdf(str(output_path))

        assert str(output_path) in written_files
        assert str(ancillary_path) in written_files

    def test_returns_none(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that NetCDF export returns None."""

        def fake_to_netcdf(self, path, *args, **kwargs):
            pass

        monkeypatch.setattr(xr.Dataset, "to_netcdf", fake_to_netcdf)

        result = l2a_product.export_to_netcdf(str(tmp_path / "scene.nc"))

        assert result is None

# ============================================================
# Test Export to Zarr
# ============================================================

class TestExportToZarr:
    
    def test_rejects_remote_paths_without_zarr_extension(
        self,
        l2a_product,
        monkeypatch,
    ):
        """Test that remote Zarr URLs require a .zarr extension."""

        monkeypatch.setattr(mod, "_zarr_major", lambda: 2)
        monkeypatch.setattr(
            l2a_product,
            "_prepare_export",
            lambda snap_compliant=False: (
                xr.Dataset(
                    {
                        "Rrs": (
                            ("y", "x"),
                            np.ones((4, 4), dtype=np.float32),
                        ),
                    },
                    coords={"x": range(4), "y": range(4)},
                ),
                {"Rrs": {"dtype": "int16"}},
            ),
        )

        with pytest.raises(ValueError, match="has to end with '.zarr'"):
            l2a_product.export_to_zarr(
                "s3://bucket/scene",
                pyramids=False,
            )


    def test_returns_store_path(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that Zarr export returns the resolved store path."""

        monkeypatch.setattr(mod, "_zarr_major", lambda: 2)
        monkeypatch.setattr(mod, "_compression_encoding", lambda *args: {})
        monkeypatch.setattr(mod, "_reset_encoding", lambda ds: None)
        monkeypatch.setattr(mod, "_load_coords", lambda ds: ds)
        monkeypatch.setattr(mod, "_clear_store", lambda *args, **kwargs: None)
        monkeypatch.setattr(mod, "_consolidate", lambda *args, **kwargs: None)

        monkeypatch.setattr(
            mod,
            "_xarray_supports_shards",
            lambda: False,
        )

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                ),
            },
            coords={"x": range(4), "y": range(4)},
            attrs={"title": "test"},
        )

        ancillary = xr.Dataset(
            {
                "aot_cams": (
                    ("yc", "xc"),
                    np.ones((2, 2), dtype=np.float32),
                ),
            },
            coords={"xc": range(2), "yc": range(2)},
        )

        monkeypatch.setattr(
            l2a_product,
            "_prepare_export",
            lambda snap_compliant=False: (
                ds.copy(),
                {"Rrs": {"dtype": "int16"}},
            ),
        )

        monkeypatch.setattr(l2a_product, "_scrub_attrs", lambda ds: None)

        l2a_product.ancillary = ancillary

        write_calls = []

        def fake_to_zarr(self, store, **kwargs):
            write_calls.append((store, kwargs))
            return MagicMock(name="delayed_write")

        monkeypatch.setattr(xr.Dataset, "to_zarr", fake_to_zarr)
        monkeypatch.setattr(
            "dask.compute",
            lambda *args, **kwargs: args,
        )

        result = l2a_product.export_to_zarr(
            str(tmp_path / "scene.zarr"),
            pyramids=False,
            sharding=False,
            consolidate=False,
            compression=None,
        )

        assert result == str(tmp_path / "scene.zarr")
        assert len(write_calls) == 3

        stores = [store for store, _ in write_calls]
        assert all(store == result for store in stores)

        groups = [kwargs.get("group") for _, kwargs in write_calls]
        assert groups == [None, "0", "ancillary"]


    def test_disables_pyramids(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that disabling pyramids writes only the full-resolution group."""

        monkeypatch.setattr(mod, "_zarr_major", lambda: 2)
        monkeypatch.setattr(mod, "_compression_encoding", lambda *args: {})
        monkeypatch.setattr(mod, "_reset_encoding", lambda ds: None)
        monkeypatch.setattr(mod, "_load_coords", lambda ds: ds)
        monkeypatch.setattr(mod, "_clear_store", lambda *args, **kwargs: None)
        monkeypatch.setattr(mod, "_xarray_supports_shards", lambda: False)

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                ),
            },
            coords={"x": range(4), "y": range(4)},
        )

        ancillary = xr.Dataset(
            {
                "aot_cams": (
                    ("yc", "xc"),
                    np.ones((2, 2), dtype=np.float32),
                ),
            },
            coords={"xc": range(2), "yc": range(2)},
        )

        monkeypatch.setattr(
            l2a_product,
            "_prepare_export",
            lambda snap_compliant=False: (
                ds.copy(),
                {"Rrs": {"dtype": "int16"}},
            ),
        )

        monkeypatch.setattr(l2a_product, "_scrub_attrs", lambda ds: None)
        l2a_product.ancillary = ancillary

        write_calls = []

        def fake_to_zarr(self, store, **kwargs):
            write_calls.append((store, kwargs))
            return MagicMock(name="delayed_write")

        monkeypatch.setattr(xr.Dataset, "to_zarr", fake_to_zarr)
        monkeypatch.setattr("dask.compute", lambda *args, **kwargs: args)

        l2a_product.export_to_zarr(
            str(tmp_path / "scene.zarr"),
            pyramids=False,
            sharding=False,
            consolidate=False,
            compression=None,
        )

        groups = [kwargs.get("group") for _, kwargs in write_calls]

        assert groups == [None, "0", "ancillary"]


    def test_forwards_storage_options(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that storage options are forwarded to Zarr writes."""

        monkeypatch.setattr(mod, "_zarr_major", lambda: 2)
        monkeypatch.setattr(mod, "_compression_encoding", lambda *args: {})
        monkeypatch.setattr(mod, "_reset_encoding", lambda ds: None)
        monkeypatch.setattr(mod, "_load_coords", lambda ds: ds)
        monkeypatch.setattr(mod, "_clear_store", lambda *args, **kwargs: None)
        monkeypatch.setattr(mod, "_xarray_supports_shards", lambda: False)

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                ),
            },
            coords={"x": range(4), "y": range(4)},
        )

        ancillary = xr.Dataset(
            {
                "aot_cams": (
                    ("yc", "xc"),
                    np.ones((2, 2), dtype=np.float32),
                ),
            },
            coords={"xc": range(2), "yc": range(2)},
        )

        monkeypatch.setattr(
            l2a_product,
            "_prepare_export",
            lambda snap_compliant=False: (
                ds.copy(),
                {"Rrs": {"dtype": "int16"}},
            ),
        )

        monkeypatch.setattr(l2a_product, "_scrub_attrs", lambda ds: None)
        l2a_product.ancillary = ancillary

        write_calls = []

        def fake_to_zarr(self, store, **kwargs):
            write_calls.append(kwargs)
            return MagicMock(name="delayed_write")

        monkeypatch.setattr(xr.Dataset, "to_zarr", fake_to_zarr)
        monkeypatch.setattr("dask.compute", lambda *args, **kwargs: args)

        storage_options = {"anon": True}

        l2a_product.export_to_zarr(
            str(tmp_path / "scene.zarr"),
            pyramids=False,
            sharding=False,
            consolidate=False,
            compression=None,
            storage_options=storage_options,
        )

        assert len(write_calls) == 3

        for kwargs in write_calls:
            assert kwargs["storage_options"] == storage_options


    def test_consolidates_metadata(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
    ):
        """Test that Zarr metadata consolidation is called after writing."""

        monkeypatch.setattr(mod, "_zarr_major", lambda: 2)
        monkeypatch.setattr(mod, "_compression_encoding", lambda *args: {})
        monkeypatch.setattr(mod, "_reset_encoding", lambda ds: None)
        monkeypatch.setattr(mod, "_load_coords", lambda ds: ds)
        monkeypatch.setattr(mod, "_clear_store", lambda *args, **kwargs: None)
        monkeypatch.setattr(mod, "_xarray_supports_shards", lambda: False)

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                ),
            },
            coords={"x": range(4), "y": range(4)},
        )

        ancillary = xr.Dataset(
            {
                "aot_cams": (
                    ("yc", "xc"),
                    np.ones((2, 2), dtype=np.float32),
                ),
            },
            coords={"xc": range(2), "yc": range(2)},
        )

        monkeypatch.setattr(
            l2a_product,
            "_prepare_export",
            lambda snap_compliant=False: (
                ds.copy(),
                {"Rrs": {"dtype": "int16"}},
            ),
        )

        monkeypatch.setattr(l2a_product, "_scrub_attrs", lambda ds: None)
        l2a_product.ancillary = ancillary

        monkeypatch.setattr(
            xr.Dataset,
            "to_zarr",
            lambda self, store, **kwargs: MagicMock(name="delayed_write"),
        )
        monkeypatch.setattr("dask.compute", lambda *args, **kwargs: args)

        mock_consolidate = MagicMock()
        monkeypatch.setattr(mod, "_consolidate", mock_consolidate)

        store = str(tmp_path / "scene.zarr")

        result = l2a_product.export_to_zarr(
            store,
            pyramids=False,
            sharding=False,
            consolidate=True,
            compression=None,
        )

        assert result == store
        mock_consolidate.assert_called_once_with(
            store,
            None,
            2,
            ["0", "ancillary"],
        )


    def test_logs_consolidation_failure(
        self,
        l2a_product,
        tmp_path,
        monkeypatch,
        caplog,
    ):
        """Test that metadata consolidation failures are logged without failing the export."""

        import logging

        monkeypatch.setattr(mod, "_zarr_major", lambda: 2)
        monkeypatch.setattr(mod, "_compression_encoding", lambda *args: {})
        monkeypatch.setattr(mod, "_reset_encoding", lambda ds: None)
        monkeypatch.setattr(mod, "_load_coords", lambda ds: ds)
        monkeypatch.setattr(mod, "_clear_store", lambda *args, **kwargs: None)
        monkeypatch.setattr(mod, "_xarray_supports_shards", lambda: False)

        ds = xr.Dataset(
            {
                "Rrs": (
                    ("y", "x"),
                    np.ones((4, 4), dtype=np.float32),
                ),
            },
            coords={"x": range(4), "y": range(4)},
        )

        ancillary = xr.Dataset(
            {
                "aot_cams": (
                    ("yc", "xc"),
                    np.ones((2, 2), dtype=np.float32),
                ),
            },
            coords={"xc": range(2), "yc": range(2)},
        )

        monkeypatch.setattr(
            l2a_product,
            "_prepare_export",
            lambda snap_compliant=False: (
                ds.copy(),
                {"Rrs": {"dtype": "int16"}},
            ),
        )

        monkeypatch.setattr(l2a_product, "_scrub_attrs", lambda ds: None)
        l2a_product.ancillary = ancillary

        monkeypatch.setattr(
            xr.Dataset,
            "to_zarr",
            lambda self, store, **kwargs: MagicMock(name="delayed_write"),
        )
        monkeypatch.setattr("dask.compute", lambda *args, **kwargs: args)

        monkeypatch.setattr(
            mod,
            "_consolidate",
            MagicMock(side_effect=RuntimeError("consolidation failed")),
        )

        with caplog.at_level(logging.WARNING):
            result = l2a_product.export_to_zarr(
                str(tmp_path / "scene.zarr"),
                pyramids=False,
                sharding=False,
                consolidate=True,
                compression=None,
            )

        assert result == str(tmp_path / "scene.zarr")
        assert "could not consolidate zarr metadata" in caplog.text
