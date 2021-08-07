r"""
These models calculate pore volumes depending on the specified shape
"""

import numpy as np
from numpy import pi as _pi


def sphere(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing pore volume values.

    """
    return _pi/6*target[pore_diameter]**3


def cube(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    Returns
    -------
    value : NumPy ndarray
        Array containing pore volume values.

    """
    return target[pore_diameter]**3


def circle(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    """
    return _pi/4 * target[pore_diameter]**2


def square(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_diameter : string
        The dictionary key of the pore diameter values

    """
    return target[pore_diameter]**2


def effective(target, pore_volume='pore.volume',
              throat_volume='throat.volume'):
    r"""
    Calculate the effective pore volume for optional use in transient
    simulations. The effective pore volume is calculated by adding half the
    volume of all neighbouring throats to the pore volume.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    pore_volume : string
        The dictionary key of the pore volume values

    throat_volume : string
        The dictionary key of the throat volume values

    Returns
    -------
    value : NumPy ndarray
        Array containing pore volume values.

    """
    network = target.project.network
    cn = network['throat.conns']
    P1 = cn[:, 0]
    P2 = cn[:, 1]
    eff_vol = np.copy(network[pore_volume])
    np.add.at(eff_vol, P1, 1/2*network[throat_volume])
    np.add.at(eff_vol, P2, 1/2*network[throat_volume])
    pores = network.map_pores(pores=target.Ps, origin=target)
    value = eff_vol[pores]
    return value
