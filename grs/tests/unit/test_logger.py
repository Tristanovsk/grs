import io
import logging
from types import SimpleNamespace
from unittest.mock import mock_open
import builtins

import pytest

from grs import class_logger as logger_module


# ============================================================
# Fixtures and helpers
# ============================================================

@pytest.fixture
def service_logger(tmp_path, monkeypatch):
    module = logger_module
    cls = module.ServiceLogger

    # Save the global logging state.
    original_methods = {
        name: getattr(logging.Logger, name)
        for name in ("error", "warning", "info", "debug")
    }

    root = logging.getLogger()
    original_handlers = root.handlers[:]
    original_level = root.level
    original_propagate = root.propagate

    original_counters = module.logcounter.copy()

    # Reset singleton state temporarily.
    monkeypatch.setattr(cls, "instance", None)
    monkeypatch.setattr(module.THIS, "klass", None)

    instance = cls(
        log_file=tmp_path / "application.log",
        error_log=tmp_path / "errors.log",
        log_level="DEBUG",
        log_console=False,
    )

    yield instance

    # Prevent close() from accessing real /proc files.
    monkeypatch.setattr(module, "get_stat", lambda *args: None)

    if getattr(instance, "first", False):
        instance.close()

    # Restore handlers.
    for handler in root.handlers[:]:
        root.removeHandler(handler)
        if handler not in original_handlers:
            handler.close()

    for handler in original_handlers:
        if handler not in root.handlers:
            root.addHandler(handler)

    root.setLevel(original_level)
    root.propagate = original_propagate

    # Restore globally modified logging methods.
    for name, method in original_methods.items():
        setattr(logging.Logger, name, method)

    module.logcounter.clear()
    module.logcounter.update(original_counters)


# ============================================================
# Tests get_instance
# ============================================================

class TestGetInstance:

    def test_returns_none_before_initialization(self,monkeypatch):
        cls = logger_module.ServiceLogger

        # Save the original state.
        old_instance = cls.instance
        old_klass = logger_module.THIS.klass

        # Temporarily reset the singleton state.
        monkeypatch.setattr(cls, "instance", None)
        monkeypatch.setattr(logger_module.THIS, "klass", None)

        assert logger_module.get_instance() is None

        # Restore the original state.
        monkeypatch.setattr(cls, "instance", old_instance)
        monkeypatch.setattr(logger_module.THIS, "klass", old_klass)



    def test_returns_logger_instance(self,service_logger):
        assert logger_module.get_instance() is service_logger


# ============================================================
# Tests Service logger configuration
# ============================================================

class TestServiceLoggerConfiguration:
    
    def test_service_logger_is_a_singleton(self,service_logger, tmp_path):
            second = logger_module.ServiceLogger(
                log_file=tmp_path / "second.log",
                error_log=tmp_path / "second_errors.log",
                log_level="INFO",
                log_console=False,
            )

            assert second is service_logger
            assert logger_module.get_instance() is service_logger
            
    def test_logger_creates_log_files(self,service_logger):
        assert service_logger.log_file is not None
        assert service_logger.error_log is not None

        assert service_logger.file_handler is not None
        assert service_logger.file_handler_error is not None

        assert service_logger.file_handler.baseFilename.endswith(
            "application.log"
        )
        assert service_logger.file_handler_error.baseFilename.endswith(
            "errors.log"
        )

    def test_logger_registers_file_handlers(self,service_logger):
        root = logging.getLogger()

        assert service_logger.file_handler in root.handlers
        assert service_logger.file_handler_error in root.handlers

    def test_logger_without_console_has_no_console_handler(self,service_logger):
        assert service_logger.console_handler is None


# ============================================================
# Tests Logging levels
# ============================================================

class TestLoggingLevels:

    def test_info_message_is_written_to_log_file(self,service_logger):
        logger = logging.getLogger("test_info")

        logger.info("Information message")

        service_logger.file_handler.flush()

        content = open(
            service_logger.log_file, encoding="utf-8"
        ).read()

        assert "INFO" in content
        assert "Information message" in content


    def test_error_message_is_written_to_both_files(self,service_logger):
        logger = logging.getLogger("test_error")

        logger.error("Something went wrong")

        service_logger.file_handler.flush()
        service_logger.file_handler_error.flush()

        general_content = open(
            service_logger.log_file, encoding="utf-8"
        ).read()
        error_content = open(
            service_logger.error_log, encoding="utf-8"
        ).read()

        assert "Something went wrong" in general_content
        assert "Something went wrong" in error_content


    def test_warning_message_is_not_written_to_error_file(self,service_logger):
        logger = logging.getLogger("test_warning")

        logger.warning("Warning message")

        service_logger.file_handler.flush()
        service_logger.file_handler_error.flush()

        error_content = open(
            service_logger.error_log, encoding="utf-8"
        ).read()

        assert "Warning message" not in error_content


    def test_log_format_contains_expected_fields(self,service_logger):
        logger = logging.getLogger("test_format")

        logger.info("Formatted message")

        service_logger.file_handler.flush()

        content = open(
            service_logger.log_file, encoding="utf-8"
        ).read()

        assert "INFO" in content
        assert "test_format" in content
        assert "Formatted message" in content
        assert "::" in content


    def test_debug_message_is_written_when_level_is_debug(self,service_logger):
        logger = logging.getLogger("test_debug")

        logger.debug("Debug message")

        service_logger.file_handler.flush()

        content = open(
            service_logger.log_file, encoding="utf-8"
        ).read()

        assert "Debug message" in content


    def test_messages_below_configured_level_are_ignored(self,tmp_path, monkeypatch):
        log_file = tmp_path / "warning.log"
        error_file = tmp_path / "warning_errors.log"

        # Ensure a fresh singleton.
        cls = logger_module.ServiceLogger
        monkeypatch.setattr(cls, "instance", None)
        monkeypatch.setattr(logger_module.THIS, "klass", None)

        logger = cls(
            log_file=log_file,
            error_log=error_file,
            log_level="WARNING",
            log_console=False,
        )

        try:
            test_logger = logging.getLogger("test_level")

            test_logger.info("Should not appear")
            test_logger.warning("Should appear")

            logger.file_handler.flush()

            content = log_file.read_text(encoding="utf-8")

            assert "Should not appear" not in content
            assert "Should appear" in content

        finally:
            # Avoid reading real /proc files during cleanup.
            monkeypatch.setattr(
                logger_module, "get_stat", lambda *args: None
            )
            logger.close()

# ============================================================
# Tests Log Counter
# ============================================================

class TestLogCounter:

    @pytest.mark.parametrize(
        ("method", "level"),
        [
            ("error", "ERROR"),
            ("warning", "WARNING"),
            ("info", "INFO"),
            ("debug", "DEBUG"),
        ],
    )
    def test_increments_for_each_level(self,service_logger, method, level):
        before = logger_module.logcounter[level]

        logger = logging.getLogger(f"test_counter_{level.lower()}")
        getattr(logger, method)("Counter test")

        assert logger_module.logcounter[level] == before + 1

    def test_does_not_increment_for_disabled_level(self,service_logger):
        logger = logging.getLogger("test_disabled")
        logger.setLevel(logging.CRITICAL + 1)

        before = logger_module.logcounter["INFO"]

        logger.info("This message is disabled")

        assert logger_module.logcounter["INFO"] == before

    def test_counters_are_independent(self,service_logger):
        before = logger_module.logcounter.copy()
        logger = logging.getLogger("test_all_counters")

        logger.error("Error")
        logger.warning("Warning")
        logger.info("Info")
        logger.debug("Debug")

        for level in ("ERROR", "WARNING", "INFO", "DEBUG"):
            assert logger_module.logcounter[level] == before[level] + 1
            
# ============================================================
# Tests Logger Output
# ============================================================

class TestLoggerOutput:

    def test_set_stat_updates_statistics(self,service_logger):
        service_logger.set_stat(
            sys_time=1.25,
            user_time=2.50,
            total_mem=1024,
        )

        assert service_logger.sys_time == 1.25
        assert service_logger.user_time == 2.50
        assert service_logger.total_mem == 1024


    def test_set_log_file_updates_filename(self,service_logger, tmp_path):
        new_path = str(tmp_path / "new.log")

        service_logger.set_log_file(new_path)

        assert service_logger.log_file == new_path


    def test_set_output_dir_updates_output_directory(self,service_logger, tmp_path):
        output_dir = str(tmp_path / "output")

        service_logger.set_output_dir(output_dir)

        assert service_logger.output_dir == output_dir


# ============================================================
# Tests Logger Close
# ============================================================

class TestLoggerClose:

    def test_clears_singleton(self,service_logger, monkeypatch):
        monkeypatch.setattr(
            logger_module, "get_stat", lambda *args: None
        )

        service_logger.close()

        # close() resets THIS.klass.
        assert logger_module.get_instance() is None

        # The class-level singleton is not reset by the implementation.
        # Reset it explicitly to avoid contaminating subsequent tests.
        monkeypatch.setattr(
            logger_module.ServiceLogger, "instance", None
        )

    def test_removes_handlers_from_root_logger(self,service_logger, monkeypatch):
        monkeypatch.setattr(
            logger_module, "get_stat", lambda *args: None
        )

        root = logging.getLogger()

        handlers = [
            service_logger.file_handler,
            service_logger.file_handler_error,
        ]

        service_logger.close()

        for handler in handlers:
            assert handler not in root.handlers


    def test_closes_file_handlers(self,service_logger, monkeypatch):
        monkeypatch.setattr(
            logger_module, "get_stat", lambda *args: None
        )

        file_handler = service_logger.file_handler
        error_handler = service_logger.file_handler_error

        service_logger.close()

        assert file_handler.stream is None
        assert error_handler.stream is None


    def test_calls_get_stat_with_current_statistics(self,service_logger, monkeypatch):
        service_logger.set_stat(1.0, 2.0, 3.0)

        calls = []
        monkeypatch.setattr(
            logger_module,
            "get_stat",
            lambda *args: calls.append(args),
        )

        service_logger.close()

        assert calls == [(1.0, 2.0, 3.0)]

    def test_propagates_get_stat_error(self,service_logger, monkeypatch):
        def fail(*args):
            raise RuntimeError("statistics unavailable")

        monkeypatch.setattr(logger_module, "get_stat", fail)

        with pytest.raises(RuntimeError, match="statistics unavailable"):
            service_logger.close()


# ============================================================
# Tests get_stat
# ============================================================

class TestGetStat:

    def test_logs_process_statistics(self,service_logger, monkeypatch):
        proc_status = "Name:\ttest\nVmPeak:\t12345 kB\n"

        stat_fields = ["0"] * 22
        stat_fields[13] = "250"
        stat_fields[14] = "150"
        stat_fields[21] = "1000"

        contents = {
            "/proc/self/status": proc_status,
            "/proc/self/stat": " ".join(stat_fields) + "\n",
            "/proc/uptime": "2000.00 1000.00\n",
        }

        original_open = builtins.open

        def fake_open(filename, *args, **kwargs):
            filename = str(filename)

            if filename in contents:
                return io.StringIO(contents[filename])

            return original_open(filename, *args, **kwargs)

        monkeypatch.setattr("builtins.open", fake_open)

        usages = [
            SimpleNamespace(
                ru_stime=1.0,
                ru_utime=2.0,
                ru_maxrss=100,
            ),
            SimpleNamespace(
                ru_stime=0.5,
                ru_utime=0.5,
                ru_maxrss=50,
            ),
        ]

        monkeypatch.setattr(
            logger_module.resource,
            "getrusage",
            lambda who: usages.pop(0),
        )

        logger_module.get_stat(
            in_sys_time=0.5,
            in_user_time=1.0,
            in_total_mem=25,
        )

        service_logger.file_handler.flush()

        content = original_open(
            service_logger.log_file, encoding="utf-8"
        ).read()

        assert "peak_vm:" in content
        assert "max_rss: 175 kB" in content
        assert "sys_cpu: 2.0 s" in content
        assert "user_cpu: 3.5 s" in content
        assert "total_run_time: 1990.0 s" in content

    def test_logs_current_counters(self,service_logger, monkeypatch):
        """Verify that get_stat includes the current logging counters."""
        import builtins
        import io

        proc_status = "VmPeak:\t1000 kB\n"

        stat_fields = ["0"] * 22
        stat_fields[13] = "100"
        stat_fields[14] = "100"
        stat_fields[21] = "100"

        contents = {
            "/proc/self/status": proc_status,
            "/proc/self/stat": " ".join(stat_fields) + "\n",
            "/proc/uptime": "200.00 100.00\n",
        }

        original_open = builtins.open

        def fake_open(filename, *args, **kwargs):
            filename = str(filename)

            if filename in contents:
                return io.StringIO(contents[filename])

            return original_open(filename, *args, **kwargs)

        monkeypatch.setattr("builtins.open", fake_open)

        usages = [
            SimpleNamespace(
                ru_stime=0.0,
                ru_utime=0.0,
                ru_maxrss=0,
            ),
            SimpleNamespace(
                ru_stime=0.0,
                ru_utime=0.0,
                ru_maxrss=0,
            ),
        ]

        monkeypatch.setattr(
            logger_module.resource,
            "getrusage",
            lambda who: usages.pop(0),
        )

        logger_module.get_stat(0, 0, 0)

        service_logger.file_handler.flush()

        with original_open(
            service_logger.log_file, encoding="utf-8"
        ) as file:
            content = file.read()

        # ERROR and WARNING counters are not modified by get_stat().
        assert f"error: {logger_module.logcounter['ERROR']}" in content
        assert f"warning: {logger_module.logcounter['WARNING']}" in content

        # get_stat() logs its statistics using logger.info(), which increments
        # the INFO counter while get_stat() is running.
        assert "info: " in content
        assert "debug: " in content