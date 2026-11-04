"""Unit tests for Kernel.final_process_chunk and kernel.final_process."""

from unittest.mock import MagicMock, Mock

import numpy as np
import pytest
import xarray as xr

from grs.grs_kernel import Kernel

# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def kernel_final_process():
    """Build a minimal kernel to test final_process_chunk."""

    kernel = Kernel.__new__(Kernel)

    # --------------------------------------------------
    # Paramètres généraux
    # --------------------------------------------------
    kernel.pressure_ref = 101320.0
    kernel.monoview = True
    kernel._cache = {}

    # Valeurs d'AOT disponibles dans les LUT
    aot_refs = np.array([0.0, 0.1, 0.2, 0.4, 0.8])
    kernel.aot_refs = aot_refs

    # Type numérique des sorties
    kernel.prod = Mock()
    kernel.prod._type = np.float32

    # --------------------------------------------------
    # Géométrie et grille
    # --------------------------------------------------
    x = np.array([100.0, 120.0])
    y = np.array([200.0, 220.0])

    wl = np.array([
        443.0, 490.0, 560.0, 665.0,
        705.0, 740.0, 783.0, 865.0, 1610.0
    ])

    # Réflectances TOA constantes
    bands = np.full((len(wl), 2, 2), 0.1, dtype=np.float32)

    angles = np.full((2, 2), 30.0, dtype=np.float32)

    raster = xr.Dataset(
        data_vars={
            "bands": (("wl", "y", "x"), bands),
            "sza": (("y", "x"), angles),
            "vza": (("y", "x"), angles),
            "raa": (("y", "x"), angles),
        },
        coords={
            "wl": wl,
            "x": x,
            "y": y,
        },
    )

    kernel.test_raster = raster

    # --------------------------------------------------
    # Carte d'AOT fournie en entrée
    # --------------------------------------------------
    kernel.test_aot_ref_raster = xr.DataArray(
        np.full((2, 2), 0.2, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
        name="aot_ref",
    )

    # --------------------------------------------------
    # Pression atmosphérique
    # --------------------------------------------------
    kernel.cams = Mock()

    kernel.cams.raster.sp = xr.DataArray(
        np.full((2, 2), 101320.0, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
    )

    # --------------------------------------------------
    # LUT et paramètres atmosphériques
    # --------------------------------------------------
    kernel.szas = np.array([30.0])
    kernel.vzas = np.array([30.0])
    kernel.azis = np.array([150.0])

    # LUT de réflectance diffuse : (wl, aot_ref, sza, vza, azi)
    rdiff_values = np.zeros(
        (len(wl), len(aot_refs), 1, 1, 1),
        dtype=np.float32,
    )

    kernel.Rdiff_lut = xr.DataArray(
        rdiff_values,
        dims=("wl", "aot_ref", "sza", "vza", "azi"),
        coords={
            "wl": wl,
            "aot_ref": aot_refs,
            "sza": [30.0],
            "vza": [30.0],
            "azi": [150.0],
        },
    )

    # LUT d'AOT spectrale
    aot_values = np.zeros(
        (len(wl), len(aot_refs)),
        dtype=np.float32,
    )

    kernel.aot_lut = xr.DataArray(
        aot_values,
        dims=("wl", "aot_ref"),
        coords={
            "wl": wl,
            "aot_ref": aot_refs,
        },
    )

    # Transmittance totale descendante
    # Les valeurs seront simulées par _interp_Tlut dans ce fixture.
    ttot_values = np.ones(
        (len(wl), len(aot_refs)),
        dtype=np.float32,
    )

    Ttot_Ed = xr.DataArray(
        ttot_values,
        dims=("wl", "aot_ref"),
        coords={
            "wl": wl,
            "aot_ref": aot_refs,
        },
    )

    kernel.trans_aero_lut = Mock()
    kernel.trans_aero_lut.Ttot_Ed = Ttot_Ed

    # --------------------------------------------------
    # Transmittances gazeuses
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
    # Rayleigh et paramètres du sunglint
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
    # Interpolation atmosphérique simulée
    # --------------------------------------------------
    kernel._R_ = Mock()

    def fake_multiplicate(spectrum, field, output):
        output[:] = spectrum[:, None, None] * field[None, :, :]
        return output

    def fake_interp_Rlut(
        szas, sza, vzas, vza, azis, azi,
        aot_refs, aot_ref_raster, Nwl, Ny, Nx, lut_values,
    ):
        return np.zeros((Nwl, Ny, Nx), dtype=np.float32)

    def fake_interp_aotlut(
        aot_ref_values, aot_ref_raster, Nwl, Ny, Nx, lut_values,
    ):
        return np.zeros((Nwl, Ny, Nx), dtype=np.float32)

    def fake_interp_Tlut(
        angles, angle_values, aot_ref_values,
        aot_ref_raster, Nwl, Ny, Nx, lut_values,
    ):
        return np.ones((Nwl, Ny, Nx), dtype=np.float32)

    kernel._R_._multiplicate.side_effect = fake_multiplicate
    kernel._R_.interp_Rlut.side_effect = fake_interp_Rlut
    kernel._R_._interp_aotlut.side_effect = fake_interp_aotlut
    kernel._R_._interp_Tlut.side_effect = fake_interp_Tlut

    # --------------------------------------------------
    # Calcul de la pente de probabilité du sunglint
    # --------------------------------------------------
    def fake_p_slope(
        sza, vza, raa, sigma2=None, monoview=True,
    ):
        return xr.DataArray(
            np.ones(sza.shape, dtype=np.float32),
            dims=sza.dims,
            coords=sza.coords,
        )

    kernel.p_slope = Mock(side_effect=fake_p_slope)

    return kernel

@pytest.fixture
def aot_ref_raster():
    return xr.DataArray(
        np.full((4, 4), 0.2, dtype=np.float32),
        dims=("y", "x"),
        coords={
            "y": np.arange(4) * 2 + 0.5,
            "x": np.arange(4) * 2 + 0.5,
        },
        name="aot_ref",
    )

@pytest.fixture
def kernel_for_final_process():
    """Kernel minimal permettant de tester final_process sans initialiser les LUT."""
    kernel = Kernel.__new__(Kernel)

    kernel.prod = MagicMock()
    kernel.prod._type = np.float32
    kernel.prod.sensor = "S2A"

    kernel.wind_img = xr.Dataset(
        {
            "wind": (
                ("y", "x"),
                np.ones((4, 4), dtype=np.float32),
            )
        },
        coords={
            "x": np.arange(4),
            "y": np.arange(4),
        },
    )

    # Le snapshot est remplacé dans les tests : aucune LUT n'est nécessaire.
    kernel._cache = {}
    return kernel


@pytest.fixture
def raster():
    """Raster minimal compatible avec la préparation de final_process."""
    wl = np.array([443, 560, 665, 865, 1610, 2200])
    y = np.arange(8)
    x = np.arange(8)

    shape = (len(wl), len(y), len(x))

    return xr.Dataset(
        {
            "bands": (
                ("wl", "y", "x"),
                np.ones(shape, dtype=np.float32),
            ),
            "sza": (
                ("y", "x"),
                np.full((8, 8), 30, dtype=np.float32),
            ),
            "raa": (
                ("y", "x"),
                np.full((8, 8), 10, dtype=np.float32),
            ),
            "vza": (
                ("y", "x"),
                np.full((8, 8), 5, dtype=np.float32),
            ),
            # Variable qui ne doit pas être transmise à map_blocks.
            "unused": (
                ("y", "x"),
                np.zeros((8, 8), dtype=np.float32),
            ),
        },
        coords={
            "wl": wl,
            "y": y,
            "x": x,
        },
    )

# ============================================================
# Tests final_process_chunk
# ============================================================

class TestFinalProcessChunk:
    
    def test_empty_block(self,kernel_final_process):
        """Returns a dataset with empty arrays when the block is empty"""
        raster = kernel_final_process.test_raster.isel(x=slice(0, 0))
        result = kernel_final_process.final_process_chunk(raster, kernel_final_process.test_aot_ref_raster)
        assert isinstance(result, xr.Dataset)
        
        assert "Rrs" in result
        assert "BRDFg" in result
        assert "aot550" in result
        
        assert result["Rrs"].shape == (9,2,0)
        assert result["BRDFg"].shape == (2,0)
        assert result["aot550"].shape == (2,0)

    def test_no_valid_pixel(self,kernel_final_process):
        """Returns empty sets when no valid pixels are available"""
        kernel = kernel_final_process
        raster = kernel_final_process.test_raster.copy(deep=True)
        raster["bands"][:] = np.nan
        
        result = kernel.final_process_chunk(raster, kernel.test_aot_ref_raster)
        
        assert isinstance(result, xr.Dataset)
            
        assert set(result.data_vars) == {"Rrs", "BRDFg", "aot550"}
        #no reflectance or sunglint has been calculated
        assert np.isnan(result["Rrs"].values).all()
        assert np.isnan(result["BRDFg"].values).all()
        #aot is interpolated independantly
        assert np.isfinite(result["aot550"].values).all()

    def test_partially_valid_pixels(self,kernel_final_process):
        """Enables compact mode when more then half of the pixels are invalid."""
        kernel = kernel_final_process
        raster = kernel.test_raster.copy(deep=True)
        
        raster["bands"][:] = np.nan
        
        #we keep only one valid pixel (y=0, x=0)
        raster["bands"].values[:,0,0] = 0.1
        
        result = kernel.final_process_chunk(raster, kernel.test_aot_ref_raster)
        
        assert isinstance(result, xr.Dataset)
        
        #the invalid pixels stay Nan for Rrs
        assert np.isnan(result["Rrs"].values[:,0,1:]).all()
        assert np.isnan(result["Rrs"].values[:,1,:]).all()
        
        #the valid pixel has been treated
        assert np.isfinite(result["Rrs"].values[:,0,0]).all()
        
        #the sunglint is calculated only for the valid pixel
        assert np.isfinite(result["BRDFg"].values[0,0])
        assert np.isnan(result["BRDFg"].values[0,1:]).all()
        assert np.isnan(result["BRDFg"].values[1, :]).all()
        
        #the aot card is finite for all pixels
        assert np.isfinite(result["aot550"].values).all()

    def test_not_compacted(self,kernel_final_process):
        """Disables compact mode when half of the pixels are valid."""
        kernel = kernel_final_process
        raster = kernel.test_raster.copy(deep=True)
        
        raster["bands"][:] = np.nan
        
        #We keep two valid pixels out of four
        raster["bands"].values[:,0,0] = 0.1    
        raster["bands"].values[:,1,0] = 0.1
            
        result = kernel.final_process_chunk(raster, kernel.test_aot_ref_raster)
        
        assert isinstance(result, xr.Dataset)
        
        #the two valid pixels were treated
        assert np.isfinite(result["Rrs"].values[:,0,0]).all()
        assert np.isfinite(result["Rrs"].values[:,1,0]).all()
        
        #the invalid pixel stay nan
        assert np.isnan(result["Rrs"].values[:,0,1]).all()
        assert np.isnan(result["Rrs"].values[:,1,1]).all()
        
        #the sunglint is calculated only for the valid pixel
        assert np.isfinite(result["BRDFg"].values[0,0]).all()
        assert np.isfinite(result["BRDFg"].values[1,0]).all()
        assert np.isnan(result["BRDFg"].values[0,1]).all()
        assert np.isnan(result["BRDFg"].values[1,1]).all()
        
        #the aot card is finite for all pixels
        assert np.isfinite(result["aot550"].values).all()
        
    def test_clips_aot_to_lut_range(self,kernel_final_process):
        """Clips the interpolated AOT to the LUT range"""
        kernel = kernel_final_process
        raster = kernel.test_raster
        # Values voluntarily outside the authorized range [0.0,0.8]    
        aot_ref_raster = xr.DataArray(
            np.array([
                [-0.5,0.2],
                [0.4,1.5]
            ]),
            dims=("y","x"),
            coords={
                "x": raster.x.values,
                "y": raster.y.values
            },
            name="aot_ref"
        )
        
        result = kernel.final_process_chunk(raster, aot_ref_raster)
        expected_aot = np.array([
            [0.0, 0.2],
            [0.4, 0.8],
        ])

        np.testing.assert_allclose(
            result["aot550"].values,
            expected_aot,
            atol=1e-6,
        )

    def test_negative_reflectance_skips_transmittance(self,kernel_final_process, monkeypatch):
        """Does not divide reflectances by the total transmittance."""
        kernel = kernel_final_process
        raster = kernel.test_raster.copy(deep=True)
        #Force negative reflectance
        raster["bands"][:] = -0.1
        
        #The total transmittance is 0.5 before condition
        monkeypatch.setattr(
            kernel._R_,
            "_interp_Tlut",
            lambda *args, **kwargs: np.full(
                (args[4], args[5], args[6]),
                0.5,
                dtype=np.float32,
            ),
        )

        result = kernel.final_process_chunk(raster, kernel.test_aot_ref_raster)
        
        #The reflectance stays negative 
        assert np.all(result["Rrs"].values < 0)

    def test_positive_reflectance_applies_transmittance(self,kernel_final_process, monkeypatch ): 
        """Divides positive reflectances by the total transmittance.""" 
        kernel = kernel_final_process 
        raster = kernel.test_raster.copy(deep=True)

        # Réflectance TOA positive et constante.
        raster["bands"][:] = 0.1

        # Simulates a null contribution from _multiplicate
        def fake_multiplicate(spectrum, field, output):
            output.fill(0)
            return output
        
        monkeypatch.setattr(
            kernel._R_,
            "_multiplicate",
            fake_multiplicate,
        )
        # On impose une transmittance totale égale à 0.5.
        monkeypatch.setattr(
            kernel._R_,
            "_interp_Tlut",
            lambda *args, **kwargs: np.full(
                (args[4], args[5], args[6]),
                0.5,
                dtype=np.float32,
            ),
        )

        result = kernel.final_process_chunk(
            raster,
            kernel.test_aot_ref_raster,
        )
        # Les réflectances finales doivent être positives.
        assert np.all(result["Rrs"].values > 0)

# ============================================================
# Tests final_process
# ============================================================

class TestFinalProcess:

    def test_creates_expected_template(self,kernel_for_final_process,raster,aot_ref_raster,monkeypatch):
        """Vérifie les dimensions, variables et type du template passé à map_blocks."""
        kernel = kernel_for_final_process

        # Le snapshot est un objet minimal dont final_process_chunk est une méthode factice.
        snapshot = MagicMock()
        snapshot.final_process_chunk = MagicMock(return_value=None)
        monkeypatch.setattr(kernel, "_snapshot", lambda: snapshot)

        captured = {}

        def fake_map_blocks(func, obj, kwargs, template):
            captured["func"] = func
            captured["obj"] = obj
            captured["kwargs"] = kwargs
            captured["template"] = template
            return template

        monkeypatch.setattr(xr, "map_blocks", fake_map_blocks)

        result = kernel.final_process(
            raster,
            aot_ref_raster,
            chunk=4,
            lazy=True,
        )

        assert set(result.data_vars) == {"Rrs", "BRDFg", "aot550"}

        assert result["Rrs"].dims == ("wl", "y", "x")
        assert result["BRDFg"].dims == ("y", "x")
        assert result["aot550"].dims == ("y", "x")

        assert result.sizes["wl"] == raster.sizes["wl"]
        assert result.sizes["y"] == raster.sizes["y"]
        assert result.sizes["x"] == raster.sizes["x"]

        assert result["Rrs"].dtype == np.float32
        assert result["BRDFg"].dtype == np.float32
        assert result["aot550"].dtype == np.float32

        # Vérifie que seules les variables utiles sont transmises.
        assert set(captured["obj"].data_vars) == {
            "bands",
            "sza",
            "raa",
            "vza",
        }

        # Vérifie que le raster a été découpé en blocs de 4 x 4.
        assert captured["obj"]["bands"].chunksizes["y"] == (4, 4)
        assert captured["obj"]["bands"].chunksizes["x"] == (4, 4)
        assert captured["obj"]["bands"].chunksizes["wl"] == (raster.sizes["wl"],)

        assert captured["func"] == snapshot.final_process_chunk
        assert captured["kwargs"]["iwl_swir"] == [-2, -1]

        assert kernel.xres is result


    def test_passes_custom_parameters(self,kernel_for_final_process,raster,aot_ref_raster,monkeypatch):
        """Vérifie la transmission de chunk et iwl_swir."""
        kernel = kernel_for_final_process

        snapshot = MagicMock()
        snapshot.final_process_chunk = MagicMock(return_value=None)
        monkeypatch.setattr(kernel, "_snapshot", lambda: snapshot)

        captured = {}

        def fake_map_blocks(func, obj, kwargs, template):
            captured["obj"] = obj
            captured["kwargs"] = kwargs
            return template

        monkeypatch.setattr(xr, "map_blocks", fake_map_blocks)

        result = kernel.final_process(
            raster,
            aot_ref_raster,
            iwl_swir=[-3, -1],
            chunk=2,
            lazy=True,
        )

        assert captured["kwargs"]["iwl_swir"] == [-3, -1]
        assert captured["obj"]["bands"].chunksizes["y"] == (2, 2, 2, 2)
        assert captured["obj"]["bands"].chunksizes["x"] == (2, 2, 2, 2)

        assert result is kernel.xres


    def test_computes_result_when_lazy_false(self,kernel_for_final_process,raster,aot_ref_raster,monkeypatch):
        """Vérifie que lazy=False déclenche compute sur le résultat."""
        kernel = kernel_for_final_process

        template = xr.Dataset(
            {
                "Rrs": (
                    ("wl", "y", "x"),
                    np.ones((6, 8, 8), dtype=np.float32),
                ),
                "BRDFg": (
                    ("y", "x"),
                    np.ones((8, 8), dtype=np.float32),
                ),
                "aot550": (
                    ("y", "x"),
                    np.ones((8, 8), dtype=np.float32),
                ),
            },
            coords={
                "wl": raster.wl,
                "y": raster.y,
                "x": raster.x,
            },
        )

        snapshot = MagicMock()
        snapshot.final_process_chunk = MagicMock(return_value=None)
        monkeypatch.setattr(kernel, "_snapshot", lambda: snapshot)

        result_mock = MagicMock()
        result_mock.compute.return_value = template

        def fake_map_blocks(*args, **kwargs):
            return result_mock

        monkeypatch.setattr(xr, "map_blocks", fake_map_blocks)

        result = kernel.final_process(
            raster,
            aot_ref_raster,
            chunk=4,
            lazy=False,
        )

        result_mock.compute.assert_called_once_with()
        assert result is template
        assert kernel.xres is template


    def test_uses_snapshot(self,kernel_for_final_process,raster,aot_ref_raster,monkeypatch):
        """Vérifie que le snapshot est créé et utilisé pour les calculs par blocs."""
        kernel = kernel_for_final_process

        snapshot = MagicMock()
        snapshot.final_process_chunk = MagicMock(return_value=None)

        snapshot_factory = MagicMock(return_value=snapshot)
        monkeypatch.setattr(kernel, "_snapshot", snapshot_factory)

        def fake_map_blocks(func, obj, kwargs, template):
            assert func == snapshot.final_process_chunk
            return template

        monkeypatch.setattr(xr, "map_blocks", fake_map_blocks)

        kernel.final_process(
            raster,
            aot_ref_raster,
            chunk=4,
        )

        snapshot_factory.assert_called_once_with()