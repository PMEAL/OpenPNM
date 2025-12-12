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
from pathlib import Path

FORMAT = "%(message)s"
logging.basicConfig(
    format=FORMAT, datefmt="[%X]", handlers=[RichHandler(rich_tracebacks=True)]
)

logger = logging.getLogger("openpnm")

_pyproject = Path(__file__).parents[2] / "pyproject.toml"

if _pyproject.exists():
    with open(_pyproject, "rb") as f:
        data = _toml.load(f)
        __version__ = data["project"]["version"]
        logger.debug("Loaded version from pyproject.toml")
else:
    __version__ = _metadata.version(__package__ or __name__)
    logger.debug("Loaded version from importlib.metadata")


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

_np.seterr(divide="ignore", invalid="ignore")
