import numpy as _np
from openpnm.models.geometry import _geodocs


__all__ = ["cylinder",
           "cuboid",
           "extrusion",
           "rectangle"]


@_geodocs
def cylinder(
    network,
    throat_diameter='throat.diameter',
    throat_length='throat.length',
):
    r"""
    Calculate surface area for a cylindrical throat

    Parameters
    ----------
    %(network)s
    %(Dt)s
    %(Lt)s

    Returns
    -------
    surface_areas : ndarray
        A numpy ndarray containing throat surface area values

    """
    return _np.pi * network[throat_diameter] * network[throat_length]


@_geodocs
def cuboid(
    network,
    throat_diameter='throat.diameter',
    throat_length='throat.length',
):
    r"""
    Calculate surface area for a cuboid throat

    Parameters
    ----------
    %(network)s
    %(Dt)s
    %(Lt)s

    Returns
    -------

    """
    return 4 * network[throat_diameter] * network[throat_length]


@_geodocs
def extrusion(
    network,
    throat_perimeter='throat.perimeter',
    throat_length='throat.length',
):
    r"""
    Calculate surface area for an arbitrary shaped throat give the perimeter
    and length.

    Parameters
    ----------
    %(network)s
    %(Pt)s
    %(Lt)s

    Returns
    -------

    """
    return network[throat_perimeter] * network[throat_length]


@_geodocs
def rectangle(
    network,
    throat_length='throat.length',
):
    r"""
    Calculate surface area for a rectangular throat

    Only suitable for true 2D simulations

    Parameters
    ----------
    %(network)s
    %(Lt)s

    Returns
    -------

    """
    return 2 * network[throat_length]
