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
        ver = ver[: -len(suffix)]
    return ver


def _setup_logger_rich():
    import logging
    from rich.logging import RichHandler

    FORMAT = "%(message)s"
    logging.basicConfig(
        format=FORMAT, datefmt="[%X]", handlers=[RichHandler(rich_tracebacks=True)]
    )


_setup_logger_rich()
