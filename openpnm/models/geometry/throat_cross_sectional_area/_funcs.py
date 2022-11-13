from numpy import pi as _pi
from openpnm.models.geometry import _geodocs


__all__ = ["cylinder", "cuboid", "rectangle"]


@_geodocs
def cylinder(
    network,
    throat_diameter='throat.diameter',
):
    r"""
    Calculate throat cross-sectional area for a cylindrical throat

    Parameters
    ----------
    %(network)s
    %(Dt)s

    Returns
    -------

    """
    diams = network[throat_diameter]
    value = _pi / 4 * diams**2
    return value


@_geodocs
def cuboid(
    network,
    throat_diameter='throat.diameter',
):
    r"""
    Calculate throat cross-sectional area for a cuboid throat

    Parameters
    ----------
    %(network)s
    %(Dt)s

    Returns
    -------

    """
    diams = network[throat_diameter]
    value = (diams)**2
    return value


@_geodocs
def rectangle(
    network,
    throat_diameter='throat.diameter',
):
    r"""
    Calculate throat cross-sectional area for a rectangular throat

    Parameters
    ----------
    %(network)s
    %(Dt)s

    Returns
    -------

    """
    return network[throat_diameter]
