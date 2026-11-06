from unittest.mock import Mock

import numpy as np
import pytest
import xarray as xr
import dask.array as dsa
import builtins

from grs.grs_kernel import Kernel
from grs.grs_kernel import _find_min_cost, _aot_ref_candidates, _unique_values, _load_small

# ============================================================
# Tests _load_small
# ============================================================

class TestLoadSmall:
    
    def test_loads_small_dask_dataarray(self):
        obj = xr.DataArray(
            dsa.from_array(np.arange(10), chunks=5)
        )

        res = _load_small(obj, max_bytes=1000)

        assert not res.chunks
        np.testing.assert_array_equal(res.values, np.arange(10))

    def test_loads_small_dask_dataset(self):
        obj = xr.Dataset({
            "a": ("x", dsa.from_array(np.arange(10), chunks=5)),
            "b": ("x", dsa.from_array(np.arange(10) * 2, chunks=5)),
        })

        res = _load_small(obj, max_bytes=1000)

        assert not res.chunks
        np.testing.assert_array_equal(res["a"].values, np.arange(10))
        np.testing.assert_array_equal(res["b"].values, np.arange(10) * 2)   

    def test_keeps_large_dask_dataarray_lazy(self):
        obj = xr.DataArray(
            dsa.from_array(np.arange(100), chunks=10)
        )

        res = _load_small(obj, max_bytes=10)

        assert res is obj
        assert res.chunks is not None

    def test_returns_numpy_dataarray_unchanged(self):
        obj = xr.DataArray(np.arange(10))

        res = _load_small(obj, max_bytes=1000)

        assert res is obj


# ============================================================
# Tests _unique_values
# ============================================================

class TestUniqueValues:
    
    def test_numpy_arrays(self):
        a = np.array([3, 1, 2, 1, 3])
        b = np.array([1,20,20])
        
        res_a, res_b = _unique_values(a,b)
        np.testing.assert_array_equal(res_a, np.array([1,2,3])) 
        np.testing.assert_array_equal(res_b, np.array([1,20]))   

    def test_empty_array(self):
        res, = _unique_values(np.array([]))

        np.testing.assert_array_equal(res, np.array([]))

    def test_dask_array(self):
        a = xr.DataArray(
            dsa.from_array(
                np.array([3, 1, 2, 1, 3]),
                chunks=2,
            )
        )

        res, = _unique_values(a)

        np.testing.assert_array_equal(res, [1, 2, 3])

    def test_import_error(self,monkeypatch):
        original_import = builtins.__import__
        def mock_import(name, *args, **kwargs):
            if name == "dask" or name.startswith("dask."):
                raise ImportError("No module named 'dask'")
            return original_import(name, *args, **kwargs)
        
        monkeypatch.setattr(builtins, "__import__", mock_import)
        
        a = np.array([3, 1, 2, 1])
        b = np.array([5, 4, 5])

        result = _unique_values(a, b)
        
        np.testing.assert_array_equal(result[0], np.array([1, 2, 3]))
        np.testing.assert_array_equal(result[1], np.array([4, 5]))
    
# ============================================================
# Tests _find_min_cost
# ============================================================

class TestFindMinCost:
    
    def test_selects_minimum(self):
        values=np.array([10,20,30])
        cost = np.array([
            [[5.0,1.0]],
            [[2.0,4.0]],
            [[3.0,2.0]],
        ])
        
        best_value, best_cost = _find_min_cost.py_func(cost, values)
        
        np.testing.assert_array_equal(best_value, np.array([[20.0,10.0]]))
        np.testing.assert_array_equal(best_cost, np.array([[2.0,1.0]]))

    def test_all_NaN(self):
        values = np.array([10,20,30])
        cost = np.array([
                [[np.nan,np.nan]],
                [[np.nan,np.nan]],
                [[np.nan,np.nan]],
        ])
            
        best_value, best_cost = _find_min_cost.py_func(cost, values)
        
        assert np.isnan(best_value[0,0])
        assert np.isnan(best_cost[0,0])

    def test_ignores_nan_candidates(self):
        """NaN values are ignored"""
        values = np.array([10, 20, 30])

        cost = np.array([
            [[np.nan]],
            [[2.0]],
            [[3.0]],
        ])

        best_value, best_cost = _find_min_cost.py_func(cost, values)

        assert best_value[0] == 20
        assert not np.isnan(best_cost[0])


# ============================================================
# Tests _aot_ref_candidates
# ============================================================

class TestAotRefCandidates:
    
    def test_relative_threshold(self):
        aot_refs = np.array([1.0, 3.0, 5.0, 10.0])
        counts = np.array([0, 1, 2, 4])
        rel_thresh = 0.5

        res = _aot_ref_candidates(aot_refs, counts, rel_thresh)

        np.testing.assert_array_equal(res, np.array([1.0, 3.0]))

    def test_values_with_NaN(self):
        aot_refs = np.array([1.0,np.nan,-4.0,10.0])
        counts = np.array([0,0,2,0])
        rel_thresh = 1
        
        res = _aot_ref_candidates(aot_refs, counts, rel_thresh)
        np.testing.assert_array_equal(res, np.array([1.0, np.nan, 10.0]))

# ============================================================
# Tests _cached
# ============================================================

class TestCached:
    
    def test_computes_and_stores_value(self):
        kernel = Kernel.__new__(Kernel)
        kernel._cache = {}

        factory = Mock(return_value=42)

        result = kernel._cached("my_key", factory)

        assert result == 42
        assert kernel._cache["my_key"] == 42
        factory.assert_called_once_with()


    def test_reuses_existing_value(self):
        kernel = Kernel.__new__(Kernel)
        kernel._cache = {"my_key": 42}

        factory = Mock(return_value=99)

        result = kernel._cached("my_key", factory)

        assert result == 42
        factory.assert_not_called()
        
        
# ============================================================
# Tests _Negative_counter
# ============================================================

class TestNegativeCounter:
    
    def test_pixels_with_negative_reflectance(self):
        """Counts valid pixels for which minimal reflectance becomes negative."""

        from grs.grs_kernel import _NegativeCounter
        counter = _NegativeCounter(base_threshold=-0.0003)

        reflectance_aot_0 = np.array([
            [[0.02, 0.01]],
            [[0.03, 0.02]],
        ])

        counter.add(reflectance_aot_0)

        assert counter.counts.tolist() == [0]

        reflectance_aot_01 = np.array([
            [[-0.01, 0.01]],
            [[ 0.02, 0.02]],
        ])

        counter.add(reflectance_aot_01)

        assert counter.counts.tolist() == [0, 1]

        reflectance_aot_02 = np.array([
            [[-0.02, -0.01]],
            [[ 0.01,  0.02]],
        ])

        counter.add(reflectance_aot_02)

        assert counter.counts.tolist() == [0, 1, 2]
        
    def test_uses_only_first_seven_bands(self):
        """Ignores bands after seven first"""
        from grs.grs_kernel import _NegativeCounter
        counter = _NegativeCounter(base_threshold=-0.0003)
        
        reflectance = np.full((9,1,2),0.01)
        
        reflectance[7,0,1] = -0.5
        
        counter.add(reflectance, nwl=7)

        counter.add(reflectance)
        assert counter.counts.tolist() == [0, 1]

# ============================================================
# Tests _load_block
# ============================================================

class TestLoadBlock:
     
    def test_selects_requested_variables(self):
        raster = xr.Dataset({
            "a" : ("x", np.array([1,2,3])),
            "b" : ("x", np.array([4,5,6])),
            "c" : ("x", np.array([7,8,9])),
        })
        
        res = Kernel._load_block(raster, ["a", "c"])
        
        assert set(res.data_vars) == { "a", "c" }
        np.testing.assert_array_equal(res["a"].values, [1,2,3])
        np.testing.assert_array_equal(res["c"].values, [7,8,9])

    def test_ignores_missing_variable(self):
        raster = xr.Dataset({
            "a" : ("x", np.array([1,2,3]))
        })
        
        res = Kernel._load_block(raster, ["a", "missing"])
        
        assert set(res.data_vars) == { "a" }
        np.testing.assert_array_equal(res["a"].values, [1,2,3])

    def test_computes_dask_data(self):
        raster = xr.Dataset({
            "a": ("x", dsa.from_array(np.array([1, 2, 3]), chunks=2)),
        })

        res = Kernel._load_block(raster, ["a"])

        assert res["a"].chunks is None
        np.testing.assert_array_equal(res["a"].values, [1, 2, 3])