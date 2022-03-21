r"""
=======
OpenPNM
=======

OpenPNM is a package for performing pore network simulations of transport in
porous materials.

OpenPNM consists of several key modules. Each module is consisted of
several classes and each class is consisted of a few methods. Here, you'll
find a comprehensive documentation of the modules, classes, and finally the
methods, occasionally with basic embedded examples on how to use them.

"""

from . import utils
from . import core
from . import network
from . import geometry
from . import phase
from . import physics
from . import models
from . import algorithms
from . import solvers
from . import integrators
from . import materials
from . import topotools
from . import _skgraph
from . import io
from . import metrics
from . import contrib

from .utils import Workspace, Project

import numpy as _np
_np.seterr(divide='ignore', invalid='ignore')

__version__ = utils._get_version()
