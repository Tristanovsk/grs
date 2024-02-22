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


import sys
import os
import os.path
import unittest
import logging
from grs import class_logger
from grs import product


class TestLogger(unittest.TestCase):
    """
        class for unitary test of logger
    """

    def test_logger(self):
        """
            unitary test for test_logger
        """
        test_path = os.path.dirname(os.path.abspath(__file__))
        log_file = test_path + '/../output/log_file.log'
        odir = test_path + '/../output/'
        error_log_file = test_path + '/../output/error.log'
        # Init logger service
        class_logger.ServiceLogger(log_file=log_file, output_dir=odir, log_level='DEBUG', log_console=True)

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
        self.assertTrue(os.path.isfile(log_file))

        # Test log_file
        fichierLog = open(log_file, "r")
        # Verify content of log file
        lines = fichierLog.readlines()
        # Test message format in log file
        self.assertTrue(("INFO | TestLogger::test_logger | START of test_logger" in lines[0]))
        self.assertTrue(("DEBUG | TestLogger::test_logger | This is DEBUG log" in lines[1]))
        self.assertTrue(("INFO | TestLogger::test_logger | This is INFO log" in lines[2]))
        self.assertTrue(("WARNING | TestLogger::test_logger | This is WARNING log" in lines[3]))
        self.assertTrue(("ERROR | TestLogger::test_logger | This is ERROR log" in lines[4]))
        self.assertTrue(("INFO | TestLogger::test_logger | END of test_logger" in lines[5]))
        fichierLog.close()

        # Test error.log
        fichierLog = open(error_log_file, "r")
        # Verify content of log file
        lines = fichierLog.readlines()
        self.assertTrue(("ERROR | TestLogger::test_logger | This is ERROR log" in lines[0]))
        fichierLog.close()

        # Recall logger instance and close it
        class_logger.get_instance().close()
