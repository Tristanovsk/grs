"""Configuration: packaged defaults, user file and aerosol settings."""

import pytest

from grs.config import auto_download_enabled, load_config
from grs.grs_process import AEROSOL_SETTINGS, _aerosol_settings


def test_aerosol_defaults(isolated_config):
    settings = _aerosol_settings()
    assert set(settings) == set(AEROSOL_SETTINGS)
    assert settings['aot_resolution'] == 240
    assert settings['aot_chunk'] == 64
    assert settings['outlier_factor'] == 1.5
    assert settings['aot_smoothing_blocks'] == 3
    assert settings['fine_window'] is None


def test_arguments_override_configuration(isolated_config):
    settings = _aerosol_settings(aot_chunk=256, fine_window=7, outlier_factor=None)
    assert settings['aot_chunk'] == 256
    assert settings['fine_window'] == 7
    assert settings['outlier_factor'] == 1.5          # None: value of the configuration


def test_user_file_changes_one_setting(isolated_config):
    isolated_config.write_text('aerosol:\n  aot_resolution: 480\n')
    load_config.cache_clear()
    settings = _aerosol_settings()
    assert settings['aot_resolution'] == 480
    assert settings['aot_chunk'] == 64                # other keys of the section kept


def test_missing_aerosol_section(isolated_config):
    isolated_config.write_text('aerosol: null\n')
    load_config.cache_clear()
    with pytest.raises(KeyError):
        _aerosol_settings()


def test_auto_download_switch(isolated_config, monkeypatch):
    assert auto_download_enabled()
    monkeypatch.setenv('GRS_NO_DOWNLOAD', '1')
    assert not auto_download_enabled()
    monkeypatch.setenv('GRS_NO_DOWNLOAD', '0')
    assert auto_download_enabled()
