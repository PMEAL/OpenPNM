r"""
Utilities and helper classes/functions
======================================

This module contains two very important classes (Project and Workspace)
as well as a number of helper classes.

"""

from re import L
from ._misc import *
from ._settings import *
from ._workspace import *
from ._project import *
from ._health import *


def _get_version():
    from openpnm.__version__ import __version__ as ver
    suffix = ".dev0"
    if ver.endswith(suffix):
        ver = ver[:-len(suffix)]
    return ver


def _setup_logger():
    import os
    import logging
    # You can add info to the logger message by inserting the desired %(item)
    # For a list of available items see:
    # https://docs.python.org/3/library/logging.html#logrecord-attributes
    # NOTE: If the calling locations appears as 'root' it's because the logger
    # was not given a name in a file somewhere.  A good option is __name__.

    try:
        terminal_width = os.get_terminal_size().columns
    except OSError:
        terminal_width = 60
    column_width = terminal_width if terminal_width < 60 else 60

    log_format = \
        '-' * column_width + '\n\
        - %(levelname)-7s: %(message)s \n\
        - SOURCE : %(name)s.%(funcName)s \n\
        - TIME   : %(asctime)s\
        \n' + '-' * column_width

    logging.basicConfig(level=logging.WARNING, format=log_format)
    del log_format


def _setup_logger_rich():
    import logging
    from rich.logging import RichHandler

    FORMAT = "%(message)s"
    logging.basicConfig(
        format=FORMAT, datefmt="[%X]", handlers=[RichHandler(rich_tracebacks=True)]
    )
