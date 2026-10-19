"""Command line helpers (grs.run)."""

import os
import re

import pytest

from grs import run
from grs.run import set_ofile

from .conftest import SAFE


def test_set_ofile():
    out = set_ofile('S2A_MSIL1C_20210211T104151_N0500_R008_T31TFJ_20230608T124831.SAFE', odir='/out',
                    level_name='L2AGRS', suffix='_V3')
    assert str(out) == '/out/S2A_MSIL2AGRS_20210211T104151_N0500_R008_T31TFJ_20230608T124831_V3'


@pytest.mark.skipif(not os.path.exists(SAFE), reason='Level-1 image not on this computer')
def test_s2_cloud_cover():
    text = open(os.path.join(SAFE, 'MTD_MSIL1C.xml')).read()
    expected = float(re.search(r'<Cloud_Coverage_Assessment>\s*([\d.]+)', text).group(1)) / 100
    assert run.s2_cloud_cover(SAFE) == pytest.approx(expected)


def test_malloc_arenas_respects_environment(monkeypatch):
    calls = []

    class Libc:
        def mallopt(self, *args):
            calls.append(args)
            return 1

    monkeypatch.setattr(run.ctypes, 'CDLL', lambda *a, **k: Libc())
    monkeypatch.setenv('MALLOC_ARENA_MAX', '4')
    run.limit_malloc_arenas()
    assert calls == []
    monkeypatch.delenv('MALLOC_ARENA_MAX')
    run.limit_malloc_arenas()
    assert calls == [(-8, 2)]                          # M_ARENA_MAX = 2
