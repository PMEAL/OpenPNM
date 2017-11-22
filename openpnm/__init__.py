import sys as _sys

# Check Python version
if _sys.version_info < (3, 5):
    raise Exception('OpenPNM requires Python 3.5 or greater to run')

__version__ = '2.0.0-b'

from . import core
from . import network
from . import topotools
from . import geometry
from . import phases
from . import physics
from . import algorithms
from . import utils
from . import io
from . import materials
#from .core import Workspace
