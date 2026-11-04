import datetime
from unittest.mock import Mock
from unittest.mock import patch, MagicMock

import numpy as np
import pytest
import xarray as xr
import dask
import dask.array as dsa
import pandas as pd

from grs.product import Product, NETCDF_ENGINE

@pytest.fixture
def product_test():
    product = Product.__new__(Product)
    raster = xr.Dataset(
            data_vars={
                "bands": (
                    ("wl", "y", "x"),
                    np.arange(12, dtype=np.float32).reshape(2, 2, 3),
                ),
            },
            coords={
                "wl": [500, 600],
                "y": [10, 20],
                "x": [10, 30, 50],
            },
            attrs={
                "constellation_id": "S2",
            },
        )
    product.raster = raster
    product.y = [10, 20]
    product.x = [10, 30, 50]
    product.height = 2
    product.width = 3
    return product

def test_product_initialization_and_metadata(): 
    """Teste l'initialisation de Product, ses attributs et ses métadonnées."""

    # Arrange : création d'un Dataset minimal
    dataset = xr.Dataset(
        data_vars={
            "bands": (
                ("wl", "y", "x"),
                np.ones((3, 2, 2), dtype=np.float32),
            ),
            "vza": (
                ("y", "x"),
                np.full((2, 2), 10.0, dtype=np.float32),
            ),
            "sza": (
                ("y", "x"),
                np.full((2, 2), 30.0, dtype=np.float32),
            ),
            "wl_true": (
                ("wl",),
                np.array([490.0, 560.0, 842.0]),
            ),
            "wl_to_process": (
                ("wl_process",),
                np.array([490, 560, 842]),
            ),
        },
        coords={
            "wl": [490.5, 560.5, 842.5],
            "x": [10.0, 11.0],
            "y": [45.0, 44.0],
            "wl_process": [0, 1, 2],
        },
        attrs={
            "satellite": "S2A",
            "acquisition_date": "2024-06-15T10:30:00.000Z",
        },
    )

    # Nécessaire pour que l'accesseur .rio soit disponible
    # Act : initialisation sans dépendance à rioxarray
    rio_mock = MagicMock()

    rio_mock.transform_bounds.return_value = (
        10.0, 44.0, 12.0, 46.0
    )

    rio_mock.bounds.return_value = (
        10.0, 44.0, 12.0, 46.0
    )

    with (
        patch.object(Product, "load_auxiliary_data", autospec=True),
        patch.object(
            xr.Dataset,
            "rio",
            new_callable=lambda: property(lambda self: rio_mock),
            create=True,
        ),
    ):
        product = Product(l1c_obj=dataset)


    # Assert : attributs essentiels
    assert product.sensor == "S2A"
    assert product.date == datetime.datetime(2024, 6, 15, 10, 30, 0)

    assert product.width == 2
    assert product.height == 2

    assert product.output == "Rrs"
    assert product.rrs is True

    assert product.auxdatabase == (
        "cams-global-atmospheric-composition-forecasts"
    )
    assert product.auxdata is None
    assert product.aero_lut is None

    assert product.chunk == 512
    assert product._type is np.float32
    assert product.nodata == -999.9

    # Assert : métadonnées du Dataset
    assert product.raster.attrs["processor"] == product.processor
    assert product.raster.attrs["version"]

    processing_time = product.raster.attrs["processing_time"]
    datetime.datetime.strptime(
        processing_time,
        "%Y-%m-%dT%H:%M:%S.%f",
    )

    # Assert : coordonnées temporelles
    assert "time" in product.raster.coords

    np.testing.assert_equal(
        product.raster.coords["time"].values,
        np.datetime64("2024-06-15T10:30:00"),
    )

    # Assert : conversion des longueurs d'onde
    np.testing.assert_array_equal(
        product.raster["wl"].values,
        [490, 560, 842],
    )
    assert np.issubdtype(product.raster["wl"].dtype, np.integer)

def test_load_auxilirary_data():
    product = Product.__new__(Product)
    fake_gas_lut = MagicMock()
    fake_twv_lut = MagicMock()
    
    product.abs_gas_file = "fake_gas.nc"
    product.water_vapor_transmittance_file = "fake_water_vapor.nc"
    
    with patch("grs.product.xr.open_dataset") as mock_open_dataset:
        mock_open_dataset.side_effect = [
            fake_gas_lut,
            fake_twv_lut,
        ]
        product.load_auxiliary_data()
    
    assert product.gas_lut is fake_gas_lut
    assert product.Twv_lut is fake_twv_lut
    
    assert mock_open_dataset.call_count == 2
    
    mock_open_dataset.assert_any_call(
        "fake_gas.nc",
        engine=NETCDF_ENGINE,
    )

    mock_open_dataset.assert_any_call(
        "fake_water_vapor.nc",
        engine=NETCDF_ENGINE,
    )

def test_get_elevation(product_test):
    product = product_test
    
    product.get_elevation();
    expected = xr.DataArray(
        np.zeros((2,3), dtype=np.float64),
        coords={
            "y": [10, 20],
            "x": [10, 30, 50],
        },
        attrs={
            "description": "Digital elevation model from Copernicus30m",
            "units": "m"
        },
        name="dem",
    )
   
    assert product.elevation.shape == (2,3)
    assert product.elevation.dims == ("y", "x")
    assert product.elevation.name == "dem"
    assert product.elevation.attrs == expected.attrs
    assert isinstance(product.elevation, xr.DataArray)
    xr.testing.assert_equal(product.elevation, expected)

def test_mu_N_monoview(product_test):
    """Check mu_N for zero angles in monoview mode."""
    product = product_test
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

    result = product.mu_N(
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
    
def test_mu_N_multiview(product_test):
    """Check mu_N calculation and dimensions in multiview mode."""

    product = product_test
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

    result = product.mu_N(
        sza=sza,
        vza=vza,
        azi=azi,
        monoview=False,
    )

    # Check output dimensions.
    assert result.dims == ("wl", "y", "x")

    # With azi = 0°, cos(alpha) = cos(sza - vza).
    # Therefore, mu_N = cos((sza + vza) / 2).
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

def test_p_slope_monoview(product_test):
    product = product_test
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
    
    #cosN = np.ones((2,2))
    #thetaN = 0
    xp_slope = 1 / (np.pi * 0.02)
    expected = np.full((2,2), xp_slope, dtype=np.float64)
    
    result = product.p_slope(
        sza=sza,
        vza=vza,
        azi=azi,
        monoview=True,
    )
    np.testing.assert_array_equal(result, expected)

def test_p_slope_multiview(product_test):
    """Check p_slope calculation and dimensions in multiview mode."""

    product = product_test
    wl = [500., 600., 700., 800.]

    # Angles vary with wavelength.
    sza_values = np.full((4), 60.0)
    vza_values = np.full((4), 60.0)
    azi_values = np.full((4), 180.)
    
    shape = (4, 2, 2)

    print(sza_values)
    
    sza = xr.DataArray(
        np.full((4,2,2), 60.0),
        dims=("wl", "y", "x"),
        coords={"wl": wl},
    )
    
    vza = xr.DataArray(
        np.full((4,2,2), 60.0),
        dims=("wl", "y", "x"),
        coords={"wl": wl},
    )

    azi = xr.DataArray(
        np.full((4,2,2), 180.0),
        dims=("wl", "y", "x"),
        coords={"wl": wl},
    )

    result = product.p_slope(
        sza=sza,
        vza=vza,
        azi=azi,
        monoview=False,
    )

    # Check output dimensions.
    assert result.dims == ("wl", "y", "x")

    #cosN = np.ones((2,2))
    #thetaN = 0
    xp_slope = 1 / (np.pi * 0.02)
    expected = np.full((4,2,2), xp_slope, dtype=np.float64)

    np.testing.assert_allclose(
        result.values,
        expected,
        rtol=1e-12,
        atol=1e-12,
    )