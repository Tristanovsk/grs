from unittest.mock import Mock
from unittest.mock import patch, MagicMock

import numpy as np
import pytest
import xarray as xr
import dask
import dask.array as dsa
import pandas as pd
from types import SimpleNamespace

from grs.utils import _getnearpos, _nearest_map, _srf_supports, _band_average
from grs.utils import Rasterization, Aerosol, GaseousTransmittance

@pytest.mark.parametrize(
    "array, value, expected",
    [
        (np.array([0,10,50,100]), 9, 1),
        (np.array([0,10,50,100]), 2000, 3),
        (np.array([0,10,50,10,100]), 9, 1),
        (np.array([0,10,50,100]), np.nan, 0),
        (np.array([0,10,50,100]), 5, 0),
    ]
)
def test_getnearpos(array, value, expected):
    assert _getnearpos.py_func(array=array, value=value) == expected

def test_nearestmap():
    array = np.array(
        [
            [1, 5],
            [np.nan, 12],
        ]
    )
    axis = np.array([0,10,15])
    
    expected = np.array(
        [
            [0, 0],
            [-1, 1]
        ]
    )
    np.testing.assert_array_equal(_nearest_map.py_func(axis=axis, values=array), expected)

# ============================================================
# Tests interp Rlut
# ============================================================   

class TestInterpRlut:
    
    def test_interp_Rlut_rayleigh_mono(self):
        raster = Rasterization.__new__(Rasterization)
        szas = np.array([0.0, 10.0, 20.0], dtype=np.float32)
        vzas = np.array([0.0, 30.0, 60.0], dtype=np.float32)
        azis = np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float32)
        
        _sza = np.array([
            [ 1.0,  9.0, 19.0],
            [11.0,  np.nan, 4.0],
        ], dtype=np.float32)

        _vza = np.array([
            [ 2.0, 29.0, 55.0],
            [45.0, 10.0,  4.0],
        ], dtype=np.float32)

        _azi = np.array([
            [  5.0,  80.0, 170.0],
            [260.0, 100.0,  10.0],
        ], dtype=np.float32)
        
        Nwl = 2
        Ny, Nx = _sza.shape #(2,3)
        lut = np.empty(
            (Nwl, len(szas), len(vzas), len(azis)),
            dtype=np.float32,
        )

        for iwl in range(Nwl):
            for isza in range(len(szas)):
                for ivza in range(len(vzas)):
                    for iazi in range(len(azis)):
                        lut[iwl, isza, ivza, iazi] = (
                            1000 * iwl
                            + 100 * isza
                            + 10 * ivza
                            + iazi
                        )
        
        result = raster._interp_Rlut_rayleigh_mono.py_func(
            szas, _sza,
            vzas, _vza,
            azis, _azi,
            Nwl, Ny, Nx, lut
        )
        
        assert result.shape == (Nwl, Ny, Nx)
        assert result.dtype == np.float32
        expected = np.array([
        [
            [0.0, 111.0, 222.0],
            [113.0, np.nan, 0.0],
        ],
        [
            [1000.0, 1111.0, 1222.0],
            [1113.0, np.nan, 1000.0],
        ],
    ], dtype=np.float32)

        np.testing.assert_allclose(
            result,
            expected,
            equal_nan=True,
        )

    def test_interp_Rlut_mono(self):
        raster = Rasterization.__new__(Rasterization)
        szas = np.array([0.0, 10.0, 20.0], dtype=np.float32)
        vzas = np.array([0.0, 30.0, 60.0], dtype=np.float32)
        azis = np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float32)
        aot_refs = np.array([0.0, 0.3, 0.6], dtype=np.float32)
        
        _sza = np.array([
            [ 1.0,  9.0, 19.0],
            [11.0,  np.nan, 4.0],
        ], dtype=np.float32)

        _vza = np.array([
            [ 2.0, 29.0, 55.0],
            [45.0, 10.0,  4.0],
        ], dtype=np.float32)

        _azi = np.array([
            [  5.0,  80.0, 170.0],
            [260.0, 100.0,  10.0],
        ], dtype=np.float32)
        
        _aot_ref = np.array([
            [0.1, 0.28, 0.52],
            [0.36, 0.02, 0.75]
        ])
        
        Nwl = 2
        Ny, Nx = _sza.shape #(2,3)
        lut = np.empty(
            (len(aot_refs), Nwl, len(szas), len(vzas), len(azis)),
            dtype=np.float32,
        )

        for iaot in range(len(aot_refs)):
            for iwl in range(Nwl):
                for isza in range(len(szas)):
                    for ivza in range(len(vzas)):
                        for iazi in range(len(azis)):
                            lut[iaot, iwl, isza, ivza, iazi] = (
                                1000 * iaot 
                                + 10000 * iwl
                                + 100 * isza
                                + 10 * ivza
                                + iazi
                            )
        
        result = raster._interp_Rlut_mono.py_func(
            szas, _sza,
            vzas, _vza,
            azis, _azi,
            aot_refs, _aot_ref,
            Nwl, Ny, Nx, lut
        )
        
        assert result.shape == (Nwl, Ny, Nx)
        assert result.dtype == np.float32
        expected = np.array([
        [
            [0.0, 1111.0, 2222.0],
            [1113.0, np.nan, 2000.0],
        ],
        [
            [10000.0, 11111.0, 12222.0],
            [11113.0, np.nan, 12000.0],
        ],
    ], dtype=np.float32)

        np.testing.assert_allclose(
            result,
            expected,
            equal_nan=True,
        )
        
    def test_interp_Rlut_rayleigh(self):
        raster = Rasterization.__new__(Rasterization)
        szas = np.array([0.0, 10.0, 20.0], dtype=np.float32)
        vzas = np.array([0.0, 30.0, 60.0], dtype=np.float32)
        azis = np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float32)

        _sza = np.array([
            [1.0,  9.0, 19.0],
            [11.0,  np.nan, 4.0],
        ], dtype=np.float32)

        _vza = np.array([
            [
                [2.0, 29.0, 55.0],
                [45.0, 10.0,  4.0]
            ],
            [
                [ 38.0, np.nan, 15.0],
                [7.0, 24.0,  50.0]
            ],
        ], dtype=np.float32)

        _azi = np.array([
            [
                [5.0,  80.0, 170.0],
                [260.0, 100.0,  10.0],
            ],
            [
                [68.0, 250.0, 8.0],
                [260.0, 162.0,  21.0],
            ]
            
        ], dtype=np.float32)

        Nwl = 2
        Ny, Nx = _sza.shape #(2,3)
        lut = np.empty(
            (Nwl, len(szas), len(vzas), len(azis)),
            dtype=np.float32,
        )

        for iwl in range(Nwl):
            for isza in range(len(szas)):
                for ivza in range(len(vzas)):
                    for iazi in range(len(azis)):
                        lut[iwl, isza, ivza, iazi] = (
                            1000 * iwl
                            + 100 * isza
                            + 10 * ivza
                            + iazi
                        )

        result = raster._interp_Rlut_rayleigh.py_func(
            szas, _sza,
            vzas, _vza,
            azis, _azi,
            Nwl, Ny, Nx, lut
        )

        assert result.shape == (Nwl, Ny, Nx)
        assert result.dtype == np.float32
        expected = np.array([
            [
                [0.0, 111.0, 222.0],
                [113.0, np.nan, 0.0],
            ],
            [
                [1011.0, np.nan, 1200.0],
                [1103.0, np.nan, 1020.0],
            ],
        ], dtype=np.float32)

        np.testing.assert_allclose(
            result,
            expected,
            equal_nan=True,
        )

    def test_interp_Rlut(self):
        raster = Rasterization.__new__(Rasterization)
        szas = np.array([0.0, 10.0, 20.0], dtype=np.float32)
        vzas = np.array([0.0, 30.0, 60.0], dtype=np.float32)
        azis = np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float32)
        aot_refs = np.array([0.0, 0.3, 0.6], dtype=np.float32)
        
        _sza = np.array([
            [1.0,  9.0, 19.0],
            [11.0,  np.nan, 4.0],
        ], dtype=np.float32)

        _vza = np.array([
            [
                [2.0, 29.0, 55.0],
                [45.0, 10.0,  4.0]
            ],
            [
                [ 38.0, np.nan, 15.0],
                [7.0, 24.0,  50.0]
            ],
        ], dtype=np.float32)

        _azi = np.array([
            [
                [5.0,  80.0, 170.0],
                [260.0, 100.0,  10.0],
            ],
            [
                [68.0, 250.0, 8.0],
                [260.0, 162.0,  21.0],
            ]
            
        ], dtype=np.float32)
        
        _aot_ref = np.array([
                [0.1, 0.28, 0.52],
                [0.36, 0.02, 0.75]
        ])

        Nwl = 2
        Ny, Nx = _sza.shape #(2,3)
        lut = np.empty(
            (len(aot_refs),Nwl, len(szas), len(vzas), len(azis)),
            dtype=np.float32,
        )

        for iaot in range(len(aot_refs)):
                for iwl in range(Nwl):
                    for isza in range(len(szas)):
                        for ivza in range(len(vzas)):
                            for iazi in range(len(azis)):
                                lut[iaot, iwl, isza, ivza, iazi] = (
                                    1000 * iaot 
                                    + 10000 * iwl
                                    + 100 * isza
                                    + 10 * ivza
                                    + iazi
                                )

        result = raster._interp_Rlut.py_func(
            szas, _sza,
            vzas, _vza,
            azis, _azi,
            aot_refs, _aot_ref,
            Nwl, Ny, Nx, lut
        )

        assert result.shape == (Nwl, Ny, Nx)
        assert result.dtype == np.float32
        expected = np.array([
            [
                [0.0, 1111.0, 2222.0],
                [1113.0, np.nan, 2000.0],
            ],
            [
                [10011.0, np.nan, 12200.0],
                [11103.0, np.nan, 12020.0],
            ],
        ], dtype=np.float32)

        np.testing.assert_allclose(
            result,
            expected,
            equal_nan=True,
        )

def test_interp_Tlut():
    raster = Rasterization.__new__(Rasterization)
    szas = np.array([0.0, 10.0, 20.0], dtype=np.float32)
    aot_refs = np.array([0.0, 0.3, 0.6], dtype=np.float32)
    
    _sza = np.array([
        [1.0,  9.0, 19.0],
        [11.0,  np.nan, 4.0],
    ], dtype=np.float32)
    
    _aot_ref = np.array([
        [0.1, 0.28, 0.52],
        [0.36, 0.02, 0.75]
    ])

    Nwl = 2
    Ny, Nx = _sza.shape #(2,3)
    lut = np.empty(
        (len(aot_refs),Nwl, len(szas)),
        dtype=np.float32,
    )

    for iaot in range(len(aot_refs)):
        for iwl in range(Nwl):
            for isza in range(len(szas)):
                    lut[iaot, iwl, isza] = (
                        10* iaot 
                        + 100 * iwl
                        + isza
                    )

    result = raster._interp_Tlut.py_func(
        szas, _sza,
        aot_refs, _aot_ref,
        Nwl, Ny, Nx, lut
    )

    assert result.shape == (Nwl, Ny, Nx)
    assert result.dtype == np.float32
    expected = np.array([
        [
            [0.0, 11.0, 22.0],
            [11.0, np.nan, 20.0],
        ],
        [
            [100.0, 111.0, 122.0],
            [111.0, np.nan, 120.0],
        ],
    ], dtype=np.float32)

    np.testing.assert_allclose(
        result,
        expected,
        equal_nan=True,
    )
    
def test_interp_aotlut():
    raster = Rasterization.__new__(Rasterization)
    aot_refs = np.array([0.0, 0.3, 0.6], dtype=np.float32)
    
    _aot_ref = np.array([
        [0.1, 0.28, 0.52],
        [0.36, np.nan, 0.75]
    ])

    Nwl = 2
    Ny, Nx = _aot_ref.shape #(2,3)
    lut = np.empty(
        (len(aot_refs),Nwl),
        dtype=np.float32,
    )

    for iaot in range(len(aot_refs)):
        for iwl in range(Nwl):
            lut[iaot, iwl] = (
                    iaot 
                + 10 * iwl
            )

    result = raster._interp_aotlut.py_func(
        aot_refs, _aot_ref,
        Nwl, Ny, Nx, lut
    )

    assert result.shape == (Nwl, Ny, Nx)
    assert result.dtype == np.float32
    expected = np.array([
        [
            [0.0, 1.0, 2.0],
            [1.0, np.nan, 2.0],
        ],
        [
            [10.0, 11.0, 12.0],
            [11.0, np.nan, 12.0],
        ],
    ], dtype=np.float32)

    np.testing.assert_allclose(
        result,
        expected,
        equal_nan=True,
    )

def test_multiplicate():
    raster = Rasterization.__new__(Rasterization)
    arrwl = np.array([2, 4])
    array = np.array([[1, 2], [3, 4]])
    arresult = np.zeros((2,2,2), dtype=np.float32)
    result = raster._multiplicate.py_func(arrwl=arrwl, raster=array, arresult=arresult)
    
    expected = np.array([
        [
            [2, 4],
            [6, 8],
        ],
        [
            [4, 8],
            [12, 16],
        ],
    ], dtype=np.float32)
    
    assert result.shape == (2,2,2)
    np.testing.assert_array_equal(result, expected)

# ============================================================
# Tests Fit Aerosol
# ============================================================   

class TestAerosolFit:
    
    @pytest.mark.parametrize(
        "wl, aot, message",
        [
            (np.array([1]), np.array([1]), 'at least 3 wavelengths are needed to fit a second-order spectral AOT model'),
            (np.array([-1, 2 ,3]), np.array([1]), 'wavelengths and AOT must be strictly positive'),
            (np.array([1, 2, 3]), np.array([-1]), 'wavelengths and AOT must be strictly positive'),
        ]
    )
    def test_fit_spectral_aot_wrong_values_raises_value_error(self,wl, aot, message):
        aerosol = Aerosol.__new__(Aerosol)
        with pytest.raises(ValueError, match=message):
            aerosol.fit_spectral_aot(wl=wl, aot=aot)

    def test_fit_spectral_aot(self):
        aerosol = Aerosol.__new__(Aerosol)
        wl = np.array([400, 500, 600, 700, 800])
        a = 2.0
        b = -1.5
        c = 0.1

        # Génération des AOT à partir du modèle
        log_wl = np.log(wl)
        aot = np.exp(a + b * log_wl + c * log_wl**2)

        aerosol.fit_spectral_aot(wl=wl, aot=aot)

        expected = np.array([a, b, c])

        np.testing.assert_allclose(
            aerosol.popt,
            expected,
            rtol=1e-8,
            atol=1e-8,
        )
        
    def test_fit_aero(self):
        nCext_f = np.array([1.4, 1.2, 0.9, 0.7])
        nCext_c = np.array([0.5, 0.7, 1.0, 1.3])

        f_expected = 0.7

        naot = (
            f_expected * nCext_f
            + (1 - f_expected) * nCext_c
        )

        aerosol = Aerosol.__new__(Aerosol)

        result = aerosol.fit_aero(
            nCext_f=nCext_f,
            nCext_c=nCext_c,
            naot=naot,
        )

        np.testing.assert_allclose(
            result,
            np.array([f_expected]),
            rtol=1e-10,
            atol=1e-10,
        )

        np.testing.assert_allclose(
            aerosol.fcoef,
            np.array([f_expected]),
            rtol=1e-10,
            atol=1e-10,
        )

        assert result.shape == (1,)

# ============================================================
# Tests _srf_supports
# ============================================================   

class TestSrfSupports:
       
    def test_nominal(self):
        wl = np.array([400., 500., 600., 700., 800.])

        srf = np.array([
            [0.0, 0.5, 1.0, 0.5, 0.0],
            [np.nan, 0.2, 0.8, 0.4, np.nan],
        ])

        supports = _srf_supports(srf, wl)

        assert len(supports) == 2

        # ---------------------------------------
        # Band 1
        # ---------------------------------------
        idx0, w0 = supports[0]

        np.testing.assert_array_equal(
            idx0,
            np.array([0, 1, 2, 3, 4]),
        )

        # Coefficients trapézoïdaux :
        # [50, 100, 100, 100, 50]
        # Multipliés par la SRF :
        # [0, 50, 100, 50, 0]
        # Normalisation par la somme (200)
        expected_w0 = np.array([0., 0.25, 0.5, 0.25, 0.])

        np.testing.assert_allclose(w0, expected_w0)

        # ---------------------------------------
        # Band 2
        # ---------------------------------------
        idx1, w1 = supports[1]

        np.testing.assert_array_equal(
            idx1,
            np.array([1, 2, 3]),
        )

        # Coefficients trapézoïdaux :
        # [50, 100, 50]
        # Multipliés par la SRF :
        # [10, 80, 20]
        # Normalisation par la somme (110)
        expected_w1 = np.array([10 / 110, 80 / 110, 20 / 110])

        np.testing.assert_allclose(w1, expected_w1)

        # Les poids de chaque bande doivent sommer à 1
        np.testing.assert_allclose(np.sum(w0), 1.0)
        np.testing.assert_allclose(np.sum(w1), 1.0)

    def test_nan_values(self):
        wl = np.array([400., 500., 600., 700., 800.])

        srf = np.array([
            [np.nan, 0.5, 1.0, 0.5, np.nan],
        ])

        supports = _srf_supports(srf, wl)

        idx, weights = supports[0]

        np.testing.assert_array_equal(
            idx,
            np.array([1, 2, 3]),
        )

        # Espacement constant : 100 nm
        # Coefficients trapézoïdaux : [50, 100, 50]
        # Produits : [25, 100, 25]
        expected_weights = np.array([
            25 / 150,
            100 / 150,
            25 / 150,
        ])

        np.testing.assert_allclose(weights, expected_weights)

    def test_insufficient_samples(self):
        wl = np.array([400., 500., 600., 700.])

        srf = np.array([
            [np.nan, 1.0, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan],
        ])

        supports = _srf_supports(srf, wl)

        # Bande 1 : un seul échantillon valide
        idx0, w0 = supports[0]

        np.testing.assert_array_equal(idx0, np.array([1]))
        assert w0.shape == (1,)
        assert np.all(np.isnan(w0))

        # Bande 2 : aucun échantillon valide
        idx1, w1 = supports[1]

        assert idx1.size == 0
        assert w1.size == 0

# ============================================================
# Tests _band_average
# ============================================================   

class TestBandAverage:
    
    def test_nominal(self):
        # Coefficients d'absorption pour 4 longueurs d'onde
        k = np.array([0.1, 0.2, 0.3, 0.4])

        # Masse d'air
        air_mass = 2.0

        # Quantités de gaz à tester
        vals = np.array([0.0, 1.0, 2.0])

        # Deux bandes spectrales :
        # bande 1 : indices 0 et 1
        # bande 2 : indices 2 et 3
        supports = [
            (
                np.array([0, 1]),
                np.array([0.25, 0.75]),
            ),
            (
                np.array([2, 3]),
                np.array([0.5, 0.5]),
            ),
        ]

        # Appel de la fonction
        result = _band_average(
            k=k,
            air_mass=air_mass,
            vals=vals,
            supports=supports,
        )

        # Calcul indépendant du résultat attendu
        expected = np.empty((2, 3))

        for b, (idx, weights) in enumerate(supports):
            for j, val in enumerate(vals):
                expected[b, j] = np.sum(
                    weights * np.exp(-air_mass * k[idx] * val)
                )

        # Vérifications
        assert result.shape == (2, 3)

        np.testing.assert_allclose(
            result,
            expected,
            rtol=1e-12,
            atol=1e-12,
        )

    def test_block_processing(self):
        k = np.array([0.1, 0.2, 0.3, 0.4])
        air_mass = 2.0
        vals = np.array([0.0, 0.5, 1.0, 1.5, 2.0])

        supports = [
            (np.array([0, 1]), np.array([0.25, 0.75])),
            (np.array([2, 3]), np.array([0.5, 0.5])),
        ]

        # Calcul sans découpage pratique en petits blocs
        result_full = _band_average(
            k, air_mass, vals, supports,
            max_elements=1000,
        )

        # Calcul avec des blocs très petits
        result_blocks = _band_average(
            k, air_mass, vals, supports,
            max_elements=2,
        )

        np.testing.assert_allclose(
            result_blocks,
            result_full,
            rtol=1e-12,
            atol=1e-12,
        )

    def test_empty_support(self):
        k = np.array([0.1, 0.2, 0.3])
        air_mass = 1.0
        vals = np.array([0.0, 1.0, 2.0])

        supports = [
            (
                np.array([], dtype=int),
                np.array([], dtype=float),
            ),
        ]

        result = _band_average(
            k=k,
            air_mass=air_mass,
            vals=vals,
            supports=supports,
        )

        np.testing.assert_allclose(
            result,
            np.zeros((1, 3)),
        )

def test_srf_data_initializes_cache():
    wl_ref = np.array([400., 500., 600., 700.])
    wl_bands = np.array([450., 650.])

    srf_hr = xr.DataArray(
        np.array([
            [0.0, 1.0, 0.5, 0.0],
            [0.0, 0.5, 1.0, 0.0],
        ]),
        dims=("wl", "wl_hr"),
        coords={
            "wl": wl_bands,
            "wl_hr": wl_ref,
        },
    )

    interp_mock = Mock(return_value=srf_hr)

    # Création d'un objet sans appeler son constructeur
    obj = GaseousTransmittance.__new__(GaseousTransmittance)

    obj._srf_cache = None
    obj.gas_lut = SimpleNamespace(
        wl=SimpleNamespace(values=wl_ref)
    )
    obj.prod = SimpleNamespace(
        raster=SimpleNamespace(
            SRF=SimpleNamespace(interp=interp_mock)
        )
    )

    # Premier appel
    supports, result_wl_bands = obj._srf_data()

    # Vérification des longueurs d'onde des bandes
    np.testing.assert_array_equal(result_wl_bands, wl_bands)

    # Vérification des indices et des poids
    idx0, w0 = supports[0]
    np.testing.assert_array_equal(idx0, [0, 1, 2, 3])
    np.testing.assert_allclose(w0, [0., 2/3, 1/3, 0.])

    idx1, w1 = supports[1]
    np.testing.assert_array_equal(idx1, [0, 1, 2, 3])
    np.testing.assert_allclose(w1, [0., 1/3, 2/3, 0.])

    # Le cache doit être initialisé
    assert obj._srf_cache is not None

    # L'interpolation ne doit avoir été effectuée qu'une fois
    interp_mock.assert_called_once_with(wl_hr=wl_ref)

    # Deuxième appel : utilisation du cache
    result_cached = obj._srf_data()

    assert result_cached is obj._srf_cache
    interp_mock.assert_called_once()

# ============================================================
# Tests _unique_values
# ============================================================   

class TestUniqueValues:
    
    def test_nominal(self):
        gt = GaseousTransmittance.__new__(GaseousTransmittance)
        a = np.array([1, np.nan ,2,10,1,np.nan,3,1,2,8])
        expected= np.array([1,2,3,8,10])
        result = gt._unique_values(a)
        np.testing.assert_array_equal(result, expected)

    def test_empty_raises_value_error(self):
        gt = GaseousTransmittance.__new__(GaseousTransmittance)
        a = np.array([np.nan,np.nan])
        
        with pytest.raises(ValueError, match="no valid value to compute the transmittance for"):
            gt._unique_values(a)
        
    def test_only_one_valid(self):
        gt = GaseousTransmittance.__new__(GaseousTransmittance)
        a = np.array([1, np.nan ,np.nan])
        expected= np.array([1, 1.2])
        result = gt._unique_values(a)
        np.testing.assert_array_equal(result, expected)