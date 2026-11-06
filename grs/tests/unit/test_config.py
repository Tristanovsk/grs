import hashlib
import io
import json
from unittest.mock import MagicMock

import pytest
import yaml

from grs import config as cfg


# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def default_config(tmp_path, monkeypatch):
    """Provide an isolated packaged configuration."""
    path = tmp_path / "config.yml"
    path.write_text(
        "# Packaged GRS configuration\n\n"
        "path:\n"
        "  grsdata: null\n"
        "  toa_lut: toa.bin\n"
        "  trans_lut: trans.bin\n"
        "zenodo:\n"
        "  record: '12345'\n"
        "  url: https://zenodo.org\n"
        "  auto_download: true\n"
        "other:\n"
        "  nested:\n"
        "    value: 1\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(cfg, "default_config_file", lambda: path)
    monkeypatch.setattr(
        cfg, "default_data_root", lambda: tmp_path / "default-lut"
    )

    # The cache must not leak configuration between tests.
    cfg.load_config.cache_clear()

    yield path

    cfg.load_config.cache_clear()


@pytest.fixture
def clean_environment(monkeypatch):
    """Remove configuration-related environment variables."""
    for name in ("GRS_CONFIG", "GRS_DATA_ROOT", "XDG_CONFIG_HOME",
                 "GRS_NO_DOWNLOAD"):
        monkeypatch.delenv(name, raising=False)


@pytest.fixture
def isolated_config(default_config, clean_environment):
    """A clean configuration environment with a temporary default file."""
    return default_config


def write_user_config(path, data):
    """Write a YAML user configuration."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


def mock_http_response(payload):
    """Create an HTTP response that behaves like a real byte stream."""
    response = MagicMock()
    response.__enter__.return_value = response
    response.__exit__.return_value = False

    def read(size=-1):
        if not hasattr(read, "consumed"):
            read.consumed = False

        if read.consumed:
            return b""

        read.consumed = True
        return payload

    response.read.side_effect = read
    return response



# ============================================================
# Tests Data root and config file
# ============================================================

class TestConfigFile:
    
    def test_default_data_root_uses_platformdirs(self,monkeypatch, tmp_path):
        expected = tmp_path / "user-data"
        monkeypatch.setattr(cfg.platformdirs, "user_data_path",
                            lambda name: expected)

        assert cfg.default_data_root() == expected / "lut"


    def test_default_config_file_returns_packaged_config(self,isolated_config):
        assert cfg.default_config_file() == isolated_config


    def test_user_config_file_uses_grs_config_environment_variable(self,monkeypatch, tmp_path, clean_environment):
        custom = tmp_path / "custom.yml"
        monkeypatch.setenv("GRS_CONFIG", str(custom))

        assert cfg.user_config_file() == custom


    def test_user_config_file_expands_home_directory(self,monkeypatch, tmp_path, clean_environment):
        monkeypatch.setenv("GRS_CONFIG", "~/custom.yml")

        original_expanduser = cfg.Path.expanduser

        def fake_expanduser(path):
            if str(path).startswith("~/"):
                return tmp_path / str(path)[2:]
            return original_expanduser(path)

        monkeypatch.setattr(cfg.Path, "expanduser", fake_expanduser)

        assert cfg.user_config_file() == tmp_path / "custom.yml"

    def test_user_config_file_uses_xdg_config_home(self,monkeypatch, tmp_path, clean_environment):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

        assert cfg.user_config_file() == tmp_path / "grs" / "config.yml"


    def test_user_config_file_defaults_to_home_config_directory(self,monkeypatch, tmp_path, clean_environment):
        monkeypatch.setattr(cfg.Path, "home", lambda: tmp_path)

        assert cfg.user_config_file() == (
            tmp_path / ".config" / "grs" / "config.yml"
        )


# ============================================================
# Tests Read user config
# ============================================================

class TestReadUserConfig:
    
    def test_returns_empty_dict_when_file_is_missing(self,monkeypatch, tmp_path, clean_environment):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

        assert cfg.read_user_config() == {}


    def test_reads_yaml(self,monkeypatch, tmp_path, clean_environment):
        user_file = tmp_path / "config.yml"
        write_user_config(user_file, {"path": {"grsdata": "/custom/lut"}})
        monkeypatch.setenv("GRS_CONFIG", str(user_file))

        assert cfg.read_user_config() == {
            "path": {"grsdata": "/custom/lut"}
        }


    def test_returns_empty_dict_for_empty_file(self,monkeypatch, tmp_path, clean_environment):
        user_file = tmp_path / "config.yml"
        user_file.write_text("", encoding="utf-8")
        monkeypatch.setenv("GRS_CONFIG", str(user_file))

        assert cfg.read_user_config() == {}


    def test_raises_for_missing_explicit_file(self,monkeypatch, tmp_path, clean_environment):
        missing = tmp_path / "missing.yml"
        monkeypatch.setenv("GRS_CONFIG", str(missing))

        with pytest.raises(FileNotFoundError, match="GRS_CONFIG"):
            cfg.read_user_config()


# ============================================================
# Tests _merge
# ============================================================

class TestMerge:
    
    def test_recursively_merges_dictionaries(self):
        base = {
            "path": {"grsdata": None, "toa_lut": "toa.bin"},
            "other": {"a": 1, "b": 2},
        }
        update = {
            "path": {"grsdata": "/custom/lut"},
            "other": {"b": 3},
        }

        result = cfg._merge(base, update)

        assert result == {
            "path": {"grsdata": "/custom/lut", "toa_lut": "toa.bin"},
            "other": {"a": 1, "b": 3},
        }

        # The original dictionaries must remain unchanged.
        assert base["path"]["grsdata"] is None
        assert base["other"]["b"] == 2


    def test_replaces_non_dictionary_values(self):
        assert cfg._merge(
            {"value": {"nested": 1}},
            {"value": 42},
        ) == {"value": 42}


    def test_adds_new_keys(self):
        assert cfg._merge({"a": 1}, {"b": 2}) == {"a": 1, "b": 2}


# ============================================================
# Tests load_config
# ============================================================

class TestLoadConfig:
    
    def test_uses_packaged_defaults(self,isolated_config):
        config = cfg.load_config()

        assert config["path"]["toa_lut"] == "toa.bin"
        assert config["path"]["trans_lut"] == "trans.bin"
        assert config["path"]["grsdata"] == str(
            isolated_config.parent / "default-lut"
        )
        assert config["zenodo"]["record"] == "12345"


    def test_merges_user_configuration(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(
            user_file,
            {
                "path": {"grsdata": "/user/lut"},
                "other": {"nested": {"value": 2}},
            },
        )
        monkeypatch.setenv("GRS_CONFIG", str(user_file))
        cfg.load_config.cache_clear()

        config = cfg.load_config()

        assert config["path"]["grsdata"] == "/user/lut"
        assert config["path"]["toa_lut"] == "toa.bin"
        assert config["other"]["nested"]["value"] == 2


    def test_environment_overrides_user_file(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(user_file, {"path": {"grsdata": "/user/lut"}})
        monkeypatch.setenv("GRS_CONFIG", str(user_file))
        monkeypatch.setenv("GRS_DATA_ROOT", "/environment/lut")
        cfg.load_config.cache_clear()

        assert cfg.load_config()["path"]["grsdata"] == "/environment/lut"


    def test_expands_user_directory(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(user_file, {"path": {"grsdata": "~/my-lut"}})
        monkeypatch.setenv("GRS_CONFIG", str(user_file))
        monkeypatch.setattr(cfg.os.path, "expanduser",
                            lambda path: path.replace("~", str(tmp_path)))
        cfg.load_config.cache_clear()

        assert cfg.load_config()["path"]["grsdata"] == str(
            tmp_path / "my-lut"
        )


    def test_cache_can_be_cleared(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(user_file, {"path": {"grsdata": "/first"}})
        monkeypatch.setenv("GRS_CONFIG", str(user_file))

        first = cfg.load_config()
        assert first["path"]["grsdata"] == "/first"

        write_user_config(user_file, {"path": {"grsdata": "/second"}})

        # The cached value remains in effect until the cache is cleared.
        assert cfg.load_config()["path"]["grsdata"] == "/first"

        cfg.load_config.cache_clear()
        assert cfg.load_config()["path"]["grsdata"] == "/second"

# ============================================================
# Tests data_root_source 
# ============================================================

class TestDataRootSource:
    
    def test_reports_environment(self,isolated_config, monkeypatch):
        monkeypatch.setenv("GRS_DATA_ROOT", "/environment/lut")

        assert cfg.data_root_source() == "GRS_DATA_ROOT"

    def test_reports_user_configuration(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(user_file, {"path": {"grsdata": "/user/lut"}})
        monkeypatch.setenv("GRS_CONFIG", str(user_file))

        assert cfg.data_root_source() == "user configuration"

    def test_reports_default_folder(isolated_config):
        assert cfg.data_root_source() == "default folder"

# ============================================================
# Tests lut_path
# ============================================================

class TestLutPath:
    
    @pytest.mark.parametrize(
        ("name", "filename"),
        [
            ("toa", "toa.bin"),
            ("transmittance", "trans.bin"),
        ],
    )
    def test_returns_expected_path(self,isolated_config, name, filename):
        path = cfg.lut_path(name)

        assert path == isolated_config.parent / "default-lut" / filename

    def test_rejects_unknown_name(self,isolated_config):
        with pytest.raises(KeyError):
            cfg.lut_path("unknown")

# ============================================================
# Tests lut_status 
# ============================================================

class TestLUTStatus:

    def test_reports_missing_files(self,isolated_config):
        status = cfg.lut_status()

        assert len(status) == 2
        assert all(exists is False for _, _, exists in status)
        assert {name for name, _, _ in status} == {
            "toa", "transmittance"
        }


    def test_reports_existing_files(self,isolated_config):
        root = isolated_config.parent / "default-lut"
        root.mkdir()
        (root / "toa.bin").write_bytes(b"toa")
        (root / "trans.bin").write_bytes(b"trans")

        status = cfg.lut_status()

        assert all(exists for _, _, exists in status)

# ============================================================
# Tests check_luts
# ============================================================

class TestCheckLuts:
    
    def test_accepts_existing_files(self,tmp_path):
        first = tmp_path / "first.bin"
        second = tmp_path / "second.bin"
        first.write_bytes(b"first")
        second.write_bytes(b"second")

        cfg.check_luts([first, second])


    def test_raises_for_missing_files(self,isolated_config, tmp_path):
        missing = tmp_path / "missing.bin"

        with pytest.raises(FileNotFoundError, match="missing.bin"):
            cfg.check_luts([missing])


    def test_error_contains_download_instructions(self,isolated_config, tmp_path):
        missing = tmp_path / "missing.bin"

        with pytest.raises(FileNotFoundError) as exc:
            cfg.check_luts([missing])

        message = str(exc.value)
        assert "grs_config --download" in message
        assert "GRS_DATA_ROOT" in message

# ============================================================
# Tests auto_download_enabled
# ============================================================

class TestAutoDownloadEnabled:
    
    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            ("", True),
            ("0", True),
            ("1", False),
            ("true", False),
            ("false", False),
        ],
    )
    def test_environment_variable(self,isolated_config, monkeypatch, value, expected):
        monkeypatch.setenv("GRS_NO_DOWNLOAD", value)

        assert cfg.auto_download_enabled() is expected


    def test_respects_configuration(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(
            user_file, {"zenodo": {"auto_download": False}}
        )
        monkeypatch.setenv("GRS_CONFIG", str(user_file))
        cfg.load_config.cache_clear()

        assert cfg.auto_download_enabled() is False


    def test_environment_overrides_configuration(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(
            user_file, {"zenodo": {"auto_download": True}}
        )
        monkeypatch.setenv("GRS_CONFIG", str(user_file))
        monkeypatch.setenv("GRS_NO_DOWNLOAD", "1")
        cfg.load_config.cache_clear()

        assert cfg.auto_download_enabled() is False


    def test_defaults_to_enabled(self,isolated_config, monkeypatch):
        monkeypatch.delenv("GRS_NO_DOWNLOAD", raising=False)

        # The configuration has no zenodo section.
        isolated_config.write_text("path:\n  grsdata: null\n", encoding="utf-8")
        cfg.load_config.cache_clear()

        assert cfg.auto_download_enabled() is True


# ============================================================
# Tests ensure_luts
# ============================================================

class TestEnsureLuts:
    
    def test_does_not_download_when_all_files_exist(self,isolated_config, monkeypatch):
        paths = [isolated_config.parent / "a", isolated_config.parent / "b"]
        status = [
            ("toa", paths[0], True),
            ("transmittance", paths[1], True),
        ]

        monkeypatch.setattr(cfg, "lut_status", lambda: status)

        def unexpected_download(*args, **kwargs):
            pytest.fail("download_luts should not be called")

        monkeypatch.setattr(cfg, "download_luts", unexpected_download)

        assert cfg.ensure_luts() == paths

    def test_downloads_when_files_are_missing(self,isolated_config, monkeypatch):
        paths = [
            isolated_config.parent / "default-lut" / "toa.bin",
            isolated_config.parent / "default-lut" / "trans.bin",
        ]

        calls = []

        def fake_download():
            calls.append(True)

            # Simulate the actual download by creating the LUT files.
            for path in paths:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fake LUT content")

        monkeypatch.setattr(cfg, "auto_download_enabled", lambda: True)
        monkeypatch.setattr(cfg, "download_luts", fake_download)

        result = cfg.ensure_luts()

        assert result == paths
        assert calls == [True]
        assert all(path.is_file() for path in paths)


    def test_raises_when_download_is_disabled(self,isolated_config, monkeypatch):
        paths = [isolated_config.parent / "toa.bin",
                isolated_config.parent / "trans.bin"]

        monkeypatch.setattr(
            cfg,
            "lut_status",
            lambda: [
                ("toa", paths[0], False),
                ("transmittance", paths[1], False),
            ],
        )
        monkeypatch.setattr(cfg, "auto_download_enabled", lambda: False)

        with pytest.raises(FileNotFoundError):
            cfg.ensure_luts()


    def test_checks_files_after_download_failure(self,isolated_config, monkeypatch):
        paths = [isolated_config.parent / "toa.bin",
                isolated_config.parent / "trans.bin"]

        monkeypatch.setattr(
            cfg,
            "lut_status",
            lambda: [
                ("toa", paths[0], False),
                ("transmittance", paths[1], False),
            ],
        )
        monkeypatch.setattr(cfg, "auto_download_enabled", lambda: True)

        def failed_download():
            raise OSError("network unavailable")

        monkeypatch.setattr(cfg, "download_luts", failed_download)

        with pytest.raises(FileNotFoundError):
            cfg.ensure_luts()


# ============================================================
# Tests init_user_config
# ============================================================

class TestInitUserConfig:
    
    def test_creates_config_file(self,isolated_config, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
        data_root = tmp_path / "my-lut"

        result = cfg.init_user_config(data_root)

        assert result.is_file()

        content = result.read_text(encoding="utf-8")
        assert str(data_root.absolute()) in content
        assert "Packaged defaults" in content
        assert "grsdata:" in content


    def test_refuses_to_overwrite_existing_file(self,isolated_config, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
        user_file = cfg.user_config_file()
        user_file.parent.mkdir(parents=True)
        user_file.write_text("existing: true\n", encoding="utf-8")

        with pytest.raises(FileExistsError):
            cfg.init_user_config(tmp_path / "lut")


    def test_force_overwrites_file(self,isolated_config, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
        user_file = cfg.user_config_file()
        user_file.parent.mkdir(parents=True)
        user_file.write_text("existing: true\n", encoding="utf-8")

        result = cfg.init_user_config(tmp_path / "new-lut", force=True)

        assert result == user_file
        assert "existing: true" not in result.read_text(encoding="utf-8")
        assert str((tmp_path / "new-lut").absolute()) in result.read_text(
            encoding="utf-8"
        )


    def test_creates_default_data_folder(self,isolated_config, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
        default_root = tmp_path / "default-lut"
        monkeypatch.setattr(cfg, "default_data_root", lambda: default_root)

        cfg.init_user_config()

        assert default_root.is_dir()


    def test_invalidates_config_cache(self,isolated_config, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
        cfg.load_config()

        cfg.init_user_config(tmp_path / "new-lut")

        assert cfg.load_config()["path"]["grsdata"] == str(
            (tmp_path / "new-lut").absolute()
        )


# ============================================================
# Tests zenodo_record_url
# ============================================================

class TestZenodoRecordUrl:
    
    def test_uses_default_url(self,isolated_config):
        assert cfg.zenodo_record_url() == (
            "https://zenodo.org/api/records/12345"
        )

    def test_removes_trailing_slash(self,isolated_config, monkeypatch, tmp_path):
        user_file = tmp_path / "user.yml"
        write_user_config(
            user_file,
            {
                "zenodo": {
                    "url": "https://example.org/",
                    "record": "abc",
                }
            },
        )
        monkeypatch.setenv("GRS_CONFIG", str(user_file))
        cfg.load_config.cache_clear()

        assert cfg.zenodo_record_url() == (
            "https://example.org/api/records/abc"
        )

    def test_raises_without_record(self,isolated_config):
        isolated_config.write_text(
            "path:\n  grsdata: null\nzenodo:\n  url: https://zenodo.org\n",
            encoding="utf-8",
        )
        cfg.load_config.cache_clear()

        with pytest.raises(ValueError, match="no Zenodo record"):
            cfg.zenodo_record_url()


# ============================================================
# Tests download_luts
# ============================================================

class TestDownloadLuts:

    def test_downloads_and_verifies_files(self,isolated_config, monkeypatch, tmp_path):
        root = tmp_path / "lut"
        root.mkdir()

        paths = [root / "toa.bin", root / "trans.bin"]
        contents = [b"TOA LUT contents", b"transmittance LUT contents"]

        monkeypatch.setattr(
            cfg,
            "lut_status",
            lambda: [
                ("toa", paths[0], False),
                ("transmittance", paths[1], False),
            ],
        )
        monkeypatch.setattr(cfg, "zenodo_record_url", lambda: "https://example.org/record")

        files = []
        for path, content in zip(paths, contents):
            files.append(
                {
                    "key": path.name,
                    "checksum": "md5:" + hashlib.md5(content).hexdigest(),
                    "size": len(content),
                    "links": {"self": f"https://example.org/{path.name}"},
                }
            )

        record_response = mock_http_response(b"")
        record_response.__enter__.return_value = io.BytesIO(
            json.dumps({"files": files, "links": {}}).encode()
        )

        file_responses = [
            mock_http_response(content) for content in contents
        ]

        responses = [record_response, *file_responses]
        monkeypatch.setattr(
            cfg.urllib.request,
            "urlopen",
            lambda url: responses.pop(0),
        )

        result = cfg.download_luts(progress=False)

        assert result == paths
        assert paths[0].read_bytes() == contents[0]
        assert paths[1].read_bytes() == contents[1]
        assert not (root / "toa.bin.part").exists()
        assert not (root / "trans.bin.part").exists()

    def test_skips_existing_files_unless_forced(self,isolated_config, monkeypatch, tmp_path):
        path = tmp_path / "toa.bin"
        path.write_bytes(b"already downloaded")

        monkeypatch.setattr(
            cfg,
            "lut_status",
            lambda: [
                ("toa", path, True),
                ("transmittance", tmp_path / "trans.bin", False),
            ],
        )
        monkeypatch.setattr(cfg, "zenodo_record_url", lambda: "https://example.org/record")

        record_response = mock_http_response(b"")
        record_response.__enter__.return_value = io.BytesIO(
            json.dumps(
                {
                    "files": [],
                    "links": {},
                }
            ).encode()
        )

        monkeypatch.setattr(
            cfg.urllib.request, "urlopen", lambda url: record_response
        )

        with pytest.raises(FileNotFoundError, match="trans.bin"):
            cfg.download_luts(progress=False)

        # The existing file remains untouched.
        assert path.read_bytes() == b"already downloaded"

    def test_raises_when_file_is_absent_from_record(self,isolated_config, monkeypatch, tmp_path):
        path = tmp_path / "missing.bin"

        monkeypatch.setattr(
            cfg,
            "lut_status",
            lambda: [
                ("toa", path, False),
                ("transmittance", tmp_path / "trans.bin", False),
            ],
        )
        monkeypatch.setattr(cfg, "zenodo_record_url", lambda: "https://example.org/record")

        response = mock_http_response(b"")
        response.__enter__.return_value = io.BytesIO(
            json.dumps({"files": [], "links": {}}).encode()
        )
        monkeypatch.setattr(
            cfg.urllib.request, "urlopen", lambda url: response
        )

        with pytest.raises(FileNotFoundError, match="missing.bin"):
            cfg.download_luts(progress=False)

    def test_removes_part_file_on_checksum_mismatch(self,isolated_config, monkeypatch, tmp_path):
        path = tmp_path / "toa.bin"
        content = b"downloaded data"

        monkeypatch.setattr(
            cfg,
            "lut_status",
            lambda: [
                ("toa", path, False),
                ("transmittance", tmp_path / "trans.bin", True),
            ],
        )
        monkeypatch.setattr(cfg, "zenodo_record_url", lambda: "https://example.org/record")

        record = {
            "files": [
                {
                    "key": path.name,
                    "checksum": "md5:incorrect-checksum",
                    "size": len(content),
                    "links": {"self": "https://example.org/toa.bin"},
                }
            ],
            "links": {},
        }

        record_response = mock_http_response(b"")
        record_response.__enter__.return_value = io.BytesIO(
            json.dumps(record).encode()
        )
        file_response = mock_http_response(content)

        responses = [record_response, file_response]
        monkeypatch.setattr(
            cfg.urllib.request,
            "urlopen",
            lambda url: responses.pop(0),
        )

        with pytest.raises(IOError, match="checksum mismatch"):
            cfg.download_luts(progress=False)

        assert not path.exists()
        assert not path.with_name(path.name + ".part").exists()


# ============================================================
# Tests main
# ============================================================

class TestMain:

    def test_returns_zero_when_luts_exist(self,isolated_config, monkeypatch, tmp_path, capsys):
        root = tmp_path / "default-lut"
        root.mkdir()
        (root / "toa.bin").write_bytes(b"toa")
        (root / "trans.bin").write_bytes(b"trans")

        assert cfg.main([]) == 0

        output = capsys.readouterr().out
        assert "LUT folder:" in output
        assert "LUT toa" in output
        assert "LUT transmittance" in output


    def test_returns_one_when_luts_are_missing(self,isolated_config, capsys):
        result = cfg.main([])

        captured = capsys.readouterr()

        assert result == 1
        assert "MISSING" in captured.out
        assert "LUT file(s) missing" in captured.err


    def test_init_refuses_to_overwrite_existing_config(self,isolated_config, monkeypatch, tmp_path, capsys):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
        user_file = cfg.user_config_file()
        user_file.parent.mkdir(parents=True)
        user_file.write_text("existing: true\n", encoding="utf-8")

        result = cfg.main(["--init", str(tmp_path / "lut")])

        captured = capsys.readouterr()

        assert result == 1
        assert "already exists" in captured.err


    def test_download_returns_one_when_download_fails(self,isolated_config, monkeypatch, capsys):
        def failed_download(force=False, progress=True):
            raise OSError("network unavailable")

        monkeypatch.setattr(cfg, "download_luts", failed_download)

        result = cfg.main(["--download"])

        captured = capsys.readouterr()

        assert result == 1
        assert "network unavailable" in captured.err