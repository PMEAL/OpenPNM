r"""
=======
OpenPNM
=======

OpenPNM is a package for performing pore network simulations of transport in
porous materials.

"""

import logging
import importlib.metadata as _metadata
import tomllib as _toml
import numpy as _np
from rich.logging import RichHandler


# try:
#     __version__ = _metadata.version(__package__ or __name__)
# except _metadata.PackageNotFoundError:
with open("./pyproject.toml", "rb") as f:
    data = _toml.load(f)
    __version__ = data["project"]["version"]

FORMAT = "%(message)s"
logging.basicConfig(
    format=FORMAT, datefmt="[%X]", handlers=[RichHandler(rich_tracebacks=True)]
)


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
