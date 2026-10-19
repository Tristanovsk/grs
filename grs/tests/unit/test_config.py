"""Packaged configuration (grs/config.yml)."""

import os

import grs
from grs.grs_process import GRSDATA, NCPU, TOALUT, TRANSLUT, Process, config


def test_configuration_keys():
    for key in ('grsdata', 'toa_lut', 'trans_lut', 'cams'):
        assert key in config['path'], key
    assert isinstance(NCPU, int) and NCPU > 0


def test_process_uses_configuration():
    process = Process()
    assert process.lut_file == os.path.join(GRSDATA, TOALUT)
    assert process.trans_lut_file == os.path.join(GRSDATA, TRANSLUT)
    assert process.Nproc == NCPU


def test_version():
    assert grs.__version__.count('.') == 2
    assert grs.__version__ in grs.__doc__            # documented in the version history
