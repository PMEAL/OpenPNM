import numpy as np
from numpy import pi as _pi
from openpnm.models.geometry import _geodocs

__all__ = ["sphere",
           "cube",
           "circle",
           "square",
           "effective"
           ]


@_geodocs
def sphere(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------
    volumes : ndarray
        Numpy ndarray containing pore volume values

    """
    return 4/3*_pi*(network[pore_diameter]/2)**3


@_geodocs
def cube(
    network,
    pore_diameter='pore.diameter'
):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    """
    return network[pore_diameter]**3


@_geodocs
def circle(
    network,
    pore_diameter='pore.diameter',
):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    """
    return _pi/4 * network[pore_diameter]**2


@_geodocs
def square(
    network,
    pore_diameter='pore.diameter',
):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    %(network)s
    %(Dp)s

    Returns
    -------

    """
    return network[pore_diameter]**2


@_geodocs
def effective(
    network,
    pore_volume='pore.volume',
    throat_volume='throat.volume',
):
    r"""
    Calculate the effective pore volume for optional use in transient
    simulations. The effective pore volume is calculated by adding half
    the volume of all neighbouring throats to the pore volume.

    Parameters
    ----------
    %(network)s
    %(Dp)s
    %(Dt)s

    Returns
    -------

    """
    cn = network['throat.conns']
    P1 = cn[:, 0]
    P2 = cn[:, 1]
    eff_vol = np.copy(network[pore_volume])
    np.add.at(eff_vol, P1, 1/2*network[throat_volume])
    np.add.at(eff_vol, P2, 1/2*network[throat_volume])
    return eff_vol
