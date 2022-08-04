r"""
Utilities and helper classes/functions
======================================

This module contains two very important classes (Project and Workspace)
as well as a number of helper classes.

"""

from ._misc import *
from ._settings import *
from ._workspace import *
from ._project import *
from ._parsers import *


def _get_version():
    from openpnm.__version__ import __version__ as ver
    suffix = ".dev0"
    if ver.endswith(suffix):
        ver = ver[:-len(suffix)]
    return ver


def _setup_logger():
    import logging
    # You can add info to the logger message by inserting the desired %(item)
    # For a list of available items see:
    # https://docs.python.org/3/library/logging.html#logrecord-attributes
    # NOTE: If the calling locations appears as 'root' it's because the logger
    # was not given a name in a file somewhere.  A good option is __name__.

    log_format = \
    '-' * 60 + '\n\
    %(levelname)-11s: %(message)s \n\
    SOURCE     : %(name)s.%(funcName)s \n\
    TIME STAMP : %(asctime)s\
    \n' + '-' * 60

    logging.basicConfig(level=logging.WARNING, format=log_format)
    del log_format
