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
"""
.. module:: class_logger.py
    :synopsis: logging class for grs
.. moduleauthor:: thales
"""

import os
import sys
import logging
import logging.handlers
import resource

# pointer to the module object instance itself.
THIS = sys.modules[__name__]

THIS.klass = None

logcounter = {'ERROR':0, 'WARNING':0, 'INFO':0, 'DEBUG':0}

# Used in unitary test
def get_instance():
    """
        This function return an instance of the class
        ServiceLogger
    """
    return THIS.klass

class ServiceLogger(logging.getLoggerClass()):
    """
        The class ServiceLogger defines all logging parameter.
        It's an interface to python logging class.
        It's define as a singleton
    """
    instance = None
    def __new__(cls, log_file: str, output_dir: str, log_level: str='INFO', log_console: bool=True):
        """
            __new__ method for class ServiceLogger
        """
        if cls.instance is None:
            cls.instance = object.__new__(cls)
        return cls.instance

    def __init__(self, log_file: str, output_dir: str, log_level: str='INFO', log_console: bool=True):
        """
            Init class ServiceLogger
        """
        THIS.klass = self
        self.output_dir = output_dir
        self.log_file = log_file
        self.sys_time = 0
        self.user_time = 0
        self.total_mem = 0

        # logging level
        # LEVEL : DEBUG, INFO, WARNING, ERROR
        # log format :
        # YYYY-MM-DDThh:mm:ss.mmm     LEVEL:ClassName:FunctionName: message
        self.log_formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d     %(levelname)s | %(name)s::%(funcName)s | %(message)s',
                                               datefmt='%Y-%m-%dT%H:%M:%S')

        # set the name of the class in log messages
        self.root_logger = logging.getLogger()
        for handler in self.root_logger.handlers:
            self.root_logger.removeHandler(handler)

        self.root_logger.propagate = False
        logging.propagate = False

        def log_error(self, message: str, *args, **kwargs):
            """
               Local function error log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            if self.isEnabledFor(logging.ERROR):
                self._log(logging.ERROR, message, args, **kwargs)
                global logcounter
                logcounter['ERROR'] += 1
        logging.Logger.error = log_error

        def log_warning(self, message: str, *args, **kwargs):
            """
               Local function warning log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            if self.isEnabledFor(logging.WARNING):
                self._log(logging.WARNING, message, args, **kwargs)
                global logcounter
                logcounter['WARNING'] += 1
        logging.Logger.warning = log_warning

        def log_info(self, message: str, *args, **kwargs):
            """
               Local function info log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            if self.isEnabledFor(logging.INFO):
                self._log(logging.INFO, message, args, **kwargs)
                global logcounter
                logcounter['INFO'] += 1
        logging.Logger.info = log_info

        def log_debug(self, message: str, *args, **kwargs):
            """
               Local function debug log

               :param message: log message
               :type message: string
               :param *args: number of argument
               :type *args: *int
               :param **kwargs: list of argument
               :type **kwargs: **type
            """
            if self.isEnabledFor(logging.DEBUG):
                self._log(logging.DEBUG, message, args, **kwargs)
                global logcounter
                logcounter['DEBUG'] += 1
        logging.Logger.debug = log_debug


        # set the logging level from the argument
        self.root_logger.setLevel(log_level)
        if not hasattr(self, 'first'):
            # First call to ServiceLogger
            self.first = True
            # create the log file
            self.file_handler = logging.FileHandler(log_file, mode='w')
            self.file_handler.setFormatter(self.log_formatter)
            self.file_handler.setLevel(log_level)
            # create the error log file
            log_error_filename = "error.log"
            log_error_file = os.path.join(output_dir, log_error_filename)
            self.file_handler_error = logging.FileHandler(log_error_file, mode='w')
            self.file_handler_error.setFormatter(self.log_formatter)
            self.file_handler_error.setLevel("ERROR")

            self.root_logger.addHandler(self.file_handler)
            self.root_logger.addHandler(self.file_handler_error)

            if log_console:
                # logging in console
                self.console_handler = logging.StreamHandler()
                self.console_handler.setFormatter(self.log_formatter)
                self.console_handler.setLevel(log_level)
                self.root_logger.addHandler(self.console_handler)
            else:
                self.console_handler = None


    def set_stat(self, sys_time:float, user_time:float, total_mem:float):
        """
           Local function debug log

           :param sys_time: sys_time of subprocess
           :type sys_time: double
           :param user_time: user_time of subprocess
           :type user_time: double
           :param total_mem: total_mem of subprocess
           :type total_mem: double
        """
        self.sys_time = sys_time
        self.user_time = user_time
        self.total_mem = total_mem

    def set_log_file(self, log_filename:str):
        """
           This method is to set the log filename

           :param log_filename: logging filename
           :type log_filename: str
        """
        self.log_file = log_filename


    def set_output_dir(self, output_dir:str):
        """
           This method is to set the output directory

           :param output_dir: output directory
           :type output_dir: str
        """
        self.output_dir = output_dir

    def close(self):
        """
            This method close the logging file
        """
        try:
            get_stat(self.sys_time, self.user_time, self.total_mem)
        except:
            raise
        self.file_handler.close()
        self.file_handler_error.close()
        self.root_logger.removeHandler(self.file_handler)
        self.root_logger.removeHandler(self.file_handler_error)
        if self.console_handler is not None:
            self.console_handler.close()
            self.root_logger.removeHandler(self.console_handler)
        del self.first
        self.instance = None
        THIS.klass = None

def get_stat(in_sys_time:float, in_user_time:float, in_total_mem:float):
    """
        This method get and log stats

        :param in_sys_time: init sys time
        :type message: float
        :param in_user_time: init user time
        :type message: float
        :param in_total_mem: init total memory
        :type message: float
    """

    # Get logger
    logger = logging.getLogger("service_logger")
    # reading /proc/self/status to get memory status
    stat_file = "/proc/self/status"
    with open(stat_file, 'r') as stat_fp:
        for ligne in stat_fp:
            if "VmPeak:" in ligne:
                logger.info(ligne.replace("VmPeak:", "peak_vm: ").replace('\t','').rstrip('\n'))
        stat_fp.close()

    # reading /proc/self/stat to get user and cpu time
    stat_file = "/proc/self/stat"
    with open(stat_file, 'r') as stat_fp:
        content = stat_fp.readlines()
        ligne = content[0].split()
        stat_fp.close()
    # Calculation of user and sys time
    user_time = float(ligne[13])/100.0
    sys_time = float(ligne[14])/100.0
    # Calculation of date start time
    process_start_time = float(ligne[21])/100.0
    # reading /proc/uptime to get system time
    stat_file = "/proc/uptime"
    with open(stat_file, 'r') as stat_fp:
        content = stat_fp.readlines()
        sys_uptime_sec = float(content[0].split()[0])
        stat_fp.close()
    # set stats in log
    resource_self = resource.getrusage(resource.RUSAGE_SELF)
    resource_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    sys_time = resource_self.ru_stime + resource_children.ru_stime + in_sys_time
    user_time = resource_self.ru_utime + resource_children.ru_utime + in_user_time
    total_mem = resource_self.ru_maxrss + resource_children.ru_maxrss + in_total_mem
    logger.info("max_rss: " + str(total_mem) + " kB")
    logger.info("sys_cpu: " + str(round(sys_time, 2)) + " s")
    logger.info("user_cpu: " + str(round(user_time, 2)) + " s")
    logger.info("total_run_time: " + str(round(float(sys_uptime_sec) - process_start_time, 2)) + " s")
    logger.info("error: " + str(logcounter['ERROR']))
    logger.info("warning: " + str(logcounter['WARNING']))
    logger.info("info: " + str(logcounter['INFO']))
    logger.info("debug: " + str(logcounter['DEBUG']))


####################################################################
####################################################################
