from numpy import pi as _pi
import numpy as _np


def sphere(target, pore_diameter='pore.diameter',
           conduit_lengths='throat.conduit_lengths'):
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
    network = target.project.network
    Rp = target[pore_diameter] / 2
    L1 = target[conduit_lengths + '.pore1']
    L2 = target[conduit_lengths + '.pore2']
    L = _np.append(L1, L2)
    am = network.create_adjacency_matrix(weights=L)
    am = _pi/3 * (-3*am.multiply(Rp[:, None]**2) + am.power(3))
    temp = _np.matlib.repeat(Rp, _np.diff(am.indptr))
    am.data += _pi*2/3 * temp**3
    V_lens = _np.array(am.sum(axis=1)).reshape(Rp.size)
    V0 = 4/3*_pi*Rp**3
    value = V0 - V_lens
    return value


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

    """
    diams = target[pore_diameter]
    value = diams**3
    return value
