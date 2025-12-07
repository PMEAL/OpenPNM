import numpy as _np
from openpnm.models.geometry import _geodocs


__all__ = ["cylinder",
           "cuboid",
           "rectangle"]


@_geodocs
def cylinder(
    network,
    throat_diameter='throat.diameter',
):
    r"""
    Calcuate the throat perimeter assuming a circular cross-section

    Parameters
    ----------
    %(network)s
    %(Dt)s

    Returns
    -------
    perimeters : ndarray
        A numpy ndarray containing throat perimeter values

    """
    return network[throat_diameter]*_np.pi


@_geodocs
def cuboid(
    network,
    throat_diameter='throat.diameter',
):
    r"""
    Calcuate the throat perimeter assuming a square cross-section

    Parameters
    ----------
    %(network)s
    %(Dt)s

    Returns
    -------

    """
    return network[throat_diameter]*4


@_geodocs
def rectangle(
    network,
    throat_diameter='throat.diameter',
):
    r"""
    Calcuate the throat perimeter assuming a rectangular cross-section (2D)

    Parameters
    ----------
    %(network)s
    %(Dt)s

    Returns
    -------

    """
    return 1.0
