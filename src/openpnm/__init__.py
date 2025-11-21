r"""
=======
OpenPNM
=======

OpenPNM is a package for performing pore network simulations of transport in
porous materials.

"""

import logging

from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    format=FORMAT, datefmt="[%X]", handlers=[RichHandler(rich_tracebacks=True)]
)

import numpy as _np

from . import (
    _skgraph,
    algorithms,
    contrib,
    core,
    integrators,
    io,
    models,
    network,
    phase,
    solvers,
    topotools,
    utils,
    visualization,
)
from .utils import Project, Workspace

_np.seterr(divide='ignore', invalid='ignore')

__version__ = utils._get_version()
