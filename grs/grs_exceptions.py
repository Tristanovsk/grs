"""
GRS specific Exceptions handler
"""


class GRSException(Exception):
    def __init__(self, message, cause=None):
        self.msg = message
        self.cause = cause

    def __str__(self):
        return repr(self.msg) + "\n" + repr(self.cause)


class GRS_IO_Exception(GRSException):
    def __init__(self, file, cause: Exception):
        super().__init__(f"I/O Exception, error while working with {file}")
        self.cause = cause
