#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# ======================================================
#
# Project : OBS2CO
#
# ======================================================
# HISTORIQUE
# FIN-HISTORIQUE
# ======================================================


import os
import os.path
import unittest
import logging
from grs import class_logger


class TestLogger(unittest.TestCase):
    """
        class for unitary test of logger
    """

    LOG_FILE = ""
    ERROR_LOG = ""

    @classmethod
    def setUpClass(cls) -> None:
        test_path = os.path.dirname(os.path.abspath(__file__))
        cls.LOG_FILE = test_path + '/../output/log_file.log'
        odir = test_path + '/../output/'
        class_logger.ServiceLogger(log_file=cls.LOG_FILE, output_dir=odir, log_level='DEBUG', log_console=True)
        cls.ERROR_LOG = test_path + '/../output/error.log'

    @classmethod
    def tearDownClass(cls) -> None:
        class_logger.get_instance().close()

    def test_logger(self):
        """
            unitary test for test_logger
        """

        # Recall the logger to get class instance
        logger = logging.getLogger(self.__class__.__name__)
        # Write some log messages
        logger.info("START of test_logger")
        # DEBUG log
        logger.debug("This is DEBUG log")
        # INFO log
        logger.info("This is INFO log")
        # WARNING log
        logger.warning("This is WARNING log")
        # ERROR log
        logger.error("This is ERROR log")
        logger.info("END of test_logger")

        # Test if log file is created
        self.assertTrue(os.path.isfile(TestLogger.LOG_FILE))

        # Test log_file
        with open(TestLogger.LOG_FILE, "r") as fichierLog:
            # Verify content of log file
            lines = fichierLog.readlines()
            # Test message format in log file
            self.assertTrue(("INFO | TestLogger::test_logger | START of test_logger" in lines[0]))
            self.assertTrue(("DEBUG | TestLogger::test_logger | This is DEBUG log" in lines[1]))
            self.assertTrue(("INFO | TestLogger::test_logger | This is INFO log" in lines[2]))
            self.assertTrue(("WARNING | TestLogger::test_logger | This is WARNING log" in lines[3]))
            self.assertTrue(("ERROR | TestLogger::test_logger | This is ERROR log" in lines[4]))
            self.assertTrue(("INFO | TestLogger::test_logger | END of test_logger" in lines[5]))

        # Test error.log
        with open(TestLogger.ERROR_LOG, "r") as fichierLog:
            # Verify content of log file
            lines = fichierLog.readlines()
            self.assertTrue(("ERROR | TestLogger::test_logger | This is ERROR log" in lines[0]))
