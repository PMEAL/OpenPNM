import numpy as _np
from openpnm.models.geometry import _geodocs
import logging


logger = logging.getLogger(__name__)


__all__ = ["sphere",
           "circle",
           "cube",
           "square"
           ]


@_geodocs
def sphere(
    network,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area',
    amin=0,
):
    r"""
    Calculates internal surface area of pore bodies assuming they are spheres,
    then subtracts the areas of the neighboring throats

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Act)s
    amin : float
        The lower limit for surface area, useful for preventing

    Returns
    -------
    surface_areas : ndarray
        Numpy ndarry containing pore surface area values

    Notes
    -----
    This function subtracts the area of the neighboring throats in a
    crude way, by simply considering the throat cross-sectional area,
    thus not accounting for the actual curvature of the intersection.

    """
    R = network[pore_diameter] / 2
    value = 4 * _np.pi * R**2
    Tca = network[throat_cross_sectional_area]
    _np.subtract.at(value, network.conns[:, 0], Tca)
    _np.subtract.at(value, network.conns[:, 1], Tca)
    if amin is not None:
        value = _np.clip(value, amin, _np.inf)
    elif _np.any(value < 0):
        logger.warn('Negative pore surface areas found')
    return value


@_geodocs
def circle(
    network,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area',
    amin=0,
):
    r"""
    Calculates internal surface area of pore bodies assuming they are
    circular then subtracts the area of the neighboring throats.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Act)s
    amin : float
        The lower limit for surface area, useful for preventing

    Returns
    -------

    Notes
    -----

    """
    value = _np.pi * network[pore_diameter]
    Tca = network[throat_cross_sectional_area]
    _np.subtract.at(value, network.conns[:, 0], Tca)
    _np.subtract.at(value, network.conns[:, 1], Tca)
    if amin is not None:
        value = _np.clip(value, amin, _np.inf)
    elif _np.any(value < 0):
        logger.warn('Negative pore surface areas found')
    return value


def cube(
    network,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area',
    amin=0,
):
    r"""
    Calculates internal surface area of pore bodies assuming they are cubes
    then subtracts the area of the neighboring throats.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Act)s
    amin : float
        The lower limit for surface area, useful for preventing

    Returns
    -------

    """
    D = network[pore_diameter]
    value = 6.0 * D**2
    Tca = network[throat_cross_sectional_area]
    _np.subtract.at(value, network.conns[:, 0], Tca)
    _np.subtract.at(value, network.conns[:, 1], Tca)
    if amin is not None:
        value = _np.clip(value, amin, _np.inf)
    elif _np.any(value < 0):
        logger.warn('Negative pore surface areas found')
    return value


def square(
    network,
    pore_diameter='pore.diameter',
    throat_cross_sectional_area='throat.cross_sectional_area',
    amin=0,
):
    r"""
    Calculates internal surface area of pore bodies assuming they are
    squares then subtracts the area of the neighboring throats.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Act)s
    amin : float
        The lower limit for surface area, useful for preventing

    Returns
    -------

    """
    D = network[pore_diameter]
    value = 4.0 * D
    Tca = network[throat_cross_sectional_area]
    _np.subtract.at(value, network.conns[:, 0], Tca)
    _np.subtract.at(value, network.conns[:, 1], Tca)
    if amin is not None:
        value = _np.clip(value, amin, _np.inf)
    elif _np.any(value < 0):
        logger.warn('Negative pore surface areas found')
    return value
