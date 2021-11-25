r"""
=======
OpenPNM
=======

OpenPNM is a package for performing pore network simulations of transport in
porous materials.

"""

from .__version__ import __version__

from . import core
from . import utils
from . import network
from . import geometry
from . import phases
from . import physics
from . import models
from . import solvers
from . import integrators
from . import algorithms
from . import materials
from . import topotools
from . import io
from . import metrics

from .utils import Workspace, Project

import numpy as _np
_np.seterr(divide='ignore', invalid='ignore')
