import sys as _sys
import logging as _logging

# Check Python version
if _sys.version_info < (3, 5):
    raise Exception('OpenPNM requires Python 3.5 or greater to run')

__version__ = '2.0.0-b'

from . import utils
from .utils import Workspace, Project
from . import core
from . import models
from . import network
from . import topotools
from . import geometry
from . import phases
from . import physics
from . import algorithms
from . import io
from . import materials


# Set up logging to file - see previous section for more details
log_format = \
    '%(asctime)s | %(levelname)-8s | %(name)s.%(funcName)s | %(message)s'
_logging.basicConfig(level=_logging.WARNING, format=log_format)
del log_format
