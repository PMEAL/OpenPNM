import numpy as np
from numpy import pi as _pi
from openpnm.utils import Docorator


docstr = Docorator()

@docstr.get_sections(base='models.geometry.pore_volume',
                     sections=['Parameters', 'Returns'])
def sphere(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    target : OpenPNM Base object
        Object with which this model is associated. This controls
        the length of the calculated array, and also provides access to
        other necessary properties.
    pore_diameter : str
        Name of the dictionary key on ``target`` where the array containing
        pore diameter values is stored

    Returns
    -------
    volumes : ndarray
        Numpy ndarray containing pore volume values

    """
    return _pi/6*target[pore_diameter]**3


@docstr.dedent
def cube(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    %(models.geometry.pore_volume.parameters)s

    Returns
    -------
    %(models.geometry.pore_volume.returns)s

    """
    return target[pore_diameter]**3


def circle(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a spherical pore body

    Parameters
    ----------
    %(models.geometry.pore_volume.parameters)s

    Returns
    -------
    %(models.geometry.pore_volume.returns)s

    """
    return _pi/4 * target[pore_diameter]**2


def square(target, pore_diameter='pore.diameter'):
    r"""
    Calculate pore volume from diameter assuming a cubic pore body

    Parameters
    ----------
    %(models.geometry.pore_volume.parameters)s

    Returns
    -------
    %(models.geometry.pore_volume.returns)s

    """
    return target[pore_diameter]**2


def effective(target, pore_volume='pore.volume',
              throat_volume='throat.volume'):
    r"""
    Calculate the effective pore volume for optional use in transient
    simulations. The effective pore volume is calculated by adding half
    the volume of all neighbouring throats to the pore volume.

    Parameters
    ----------
    %(models.geometry.pore_volume.parameters)s
    pore_volume : str
        Name of the dictionary key on ``target`` where the array containing
        pore volume values is stored
    throat_volume : string
        Name of the dictionary key on ``target`` where the array containing
        throat volume values is stored

    Returns
    -------
    %(models.geometry.pore_volume.returns)s

    """
    network = target.project.network
    cn = network['throat.conns']
    P1 = cn[:, 0]
    P2 = cn[:, 1]
    eff_vol = np.copy(network[pore_volume])
    np.add.at(eff_vol, P1, 1/2*network[throat_volume])
    np.add.at(eff_vol, P2, 1/2*network[throat_volume])
    pores = network.pores(target.name)
    value = eff_vol[pores]
    return value
