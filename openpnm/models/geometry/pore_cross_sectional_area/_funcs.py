from numpy import pi as _pi
from openpnm.models.geometry import _geodocs


__all__ = [
    "cone",
    "sphere",
    "cube",
    "circle",
    "square"
]


@_geodocs
def sphere(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate cross-sectional area assuming the pore body is a sphere

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------
    areas : ndarray
        A numpy ndarry containing pore cross-sectional area values

    """
    D = network[pore_diameter]
    return _pi/4 * D**2


@_geodocs
def cone(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate cross-sectional area assuming the pore body is a cone

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    """
    D = network[pore_diameter]
    return _pi / 4 * D**2


@_geodocs
def cube(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate cross-sectional area assuming the pore body is a cube

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    """
    D = network[pore_diameter]
    return D**2


@_geodocs
def circle(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate cross-sectional area assuming the pore body is a circle

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return network[pore_diameter]


@_geodocs
def square(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate cross-sectional area assuming the pore body is a square

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    Notes
    -----
    This model should only be used for true 2D networks, i.e. with planar
    symmetry.

    """
    return network[pore_diameter]
