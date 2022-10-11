r"""
=======
OpenPNM
=======

OpenPNM is a package for performing pore network simulations of transport in
porous materials.

"""

from . import _skgraph
from . import utils
from . import core
from . import models
from . import topotools
from . import network
from . import phase
from . import algorithms
from . import solvers
from . import integrators
from . import io
from . import contrib
from . import visualization

from .utils import Workspace, Project

import numpy as _np
_np.seterr(divide='ignore', invalid='ignore')

__version__ = utils._get_version()

utils._setup_logger_rich()
