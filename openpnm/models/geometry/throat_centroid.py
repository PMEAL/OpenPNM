import scipy as _sp
from .throat_length import ctc


def straight(target, pore_diameter='pore.diameter',
             throat_diameter='throat.diameter',
             throat_length='throat.length'):
    r"""
    Calculate coordinates of throat centroid, assuming throats don't overlap
    with their adjacent pores.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    throat_length : string
        Dictionary key of the throat length values

    Notes
    -----
    This model should only be used for Cubic networks without diagonal
    connections.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = ctc(target=target, pore_diameter=pore_diameter)
    unit_vec = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L
    Lt = network[throat_length][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    return xyz[cn[:, 0]] + (D1+Lt)/2 * unit_vec


def spherical_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter',
                    throat_length='throat.length'):
    r"""
    Calculate the coordinates of throat centroid, assuming spherical pores.
    This model accounts for the overlapping lens between pores and throats.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    throat_diameter : string
        Dictionary key of the throat diameter values

    throat_length : string
        Dictionary key of the throat length values

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = ctc(target=target, pore_diameter=pore_diameter)
    unit_vec = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L.reshape((len(L), 1))
    Lt = network[throat_length][throats]
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    L1 = _sp.zeros_like(D1)
    # Handle cases where throat diam > pore diam
    mask = Dt > D1
    L1[mask] = D1[mask]/2
    # Calcuate effective length of pore 1
    L1[~mask] = _sp.sqrt(D1[~mask]**2 - Dt[~mask]**2)/2
    midpoint = xyz[cn[:, 0]] + (L1 + Lt/2).reshape((len(L1), 1)) * unit_vec
    return midpoint
