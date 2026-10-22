"""Log files of the processing (grs.class_logger.ServiceLogger)."""

import logging

import pytest

from grs import class_logger


@pytest.fixture
def service_logger(tmp_path):
    root = logging.getLogger()
    handlers, level = root.handlers[:], root.level
    log_file, error_log = tmp_path / 'log_file.log', tmp_path / 'error.log'
    logger = class_logger.ServiceLogger(log_file=log_file, error_log=error_log, log_level='DEBUG',
                                        log_console=False)
    yield log_file, error_log
    logger.close()
    # the service logger removes the handlers of the root logger: put back those of pytest
    for handler in handlers:
        root.addHandler(handler)
    root.setLevel(level)


def test_log_files(service_logger):
    log_file, error_log = service_logger
    logger = logging.getLogger('TestLogger')
    logger.debug('This is DEBUG log')
    logger.info('This is INFO log')
    logger.warning('This is WARNING log')
    logger.error('This is ERROR log')

    # the function name is the one of the counting wrapper of the service logger
    lines = log_file.read_text().splitlines()
    assert 'DEBUG | TestLogger::log_debug | This is DEBUG log' in lines[0]
    assert 'INFO | TestLogger::log_info | This is INFO log' in lines[1]
    assert 'WARNING | TestLogger::log_warning | This is WARNING log' in lines[2]
    assert 'ERROR | TestLogger::log_error | This is ERROR log' in lines[3]
    # only the errors in the error log
    errors = error_log.read_text().splitlines()
    assert len(errors) == 1 and 'This is ERROR log' in errors[0]
