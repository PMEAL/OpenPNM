import logging
import numpy as np
from pathlib import Path


__all__ = [
    '_parse_filename',
    '_parse_args',
]


logger = logging.getLogger(__name__)


def _parse_filename(filename, ext=""):
    p = Path(filename)
    p = p.resolve()
    if ext == "":
        ext = p.suffix
    # If extension not part of filename
    ext = "." + (ext.strip("."))
    if p.suffix != ext:
        p = p.with_suffix(ext)
    return p


def _parse_args(network, phases):
    try:
        project = network.project
        network = [network]
    except AttributeError:
        project = network[0].project
    # Ensure phases is a list, even if empty
    if not isinstance(phases, list):
        phases = [phases]
    return (project, network, phases)
