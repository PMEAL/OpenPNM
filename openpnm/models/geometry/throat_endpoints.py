import numpy as _np
from .throat_length import ctc as _ctc


def cubic_pores(target, pore_diameter='pore.diameter'):
    r"""
    Calculate coordinates of throat endpoints, assuming throats don't overlap
    with their adjacent pores. This model could be applied to conduits such as
    cuboids or cylinders in series.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    This model is only accurate for cubic networks without diagonal
    connections.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = _ctc(target=target)
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    unit_vec = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
    EP1 = xyz[cn[:, 0]] + 0.5 * D1[:, _np.newaxis] * unit_vec
    EP2 = xyz[cn[:, 1]] - 0.5 * D2[:, _np.newaxis] * unit_vec
    # Handle overlapping pores
    overlap = L - 0.5 * (D1 + D2) < 0
    mask = (D1 >= D2) & overlap
    EP2[mask] = EP1[mask]
    mask = (D1 < D2) & overlap
    EP1[mask] = EP2[mask]
    return {'head': EP1, 'tail': EP2}


def square_pores(target, pore_diameter='pore.diameter'):
    r"""
    Calculate coordinates of throat endpoints, assuming throats don't overlap
    with their adjacent pores. This model could be applied to conduits such as
    cuboids or cylinders in series in true 2D simulations.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    This model is only accurate for cubic networks without diagonal
    connections.

    """
    return cubic_pores(target, pore_diameter=pore_diameter)


def straight_throat(target, throat_centroid='throat.centroid',
                    throat_vector='throat.vector',
                    throat_length='throat.length'):
    r"""
    Calculate the coordinates of throat endpoints given a central coordinate,
    unit vector along the throat direction and a length.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_centroid : string
        Dictionary key of the throat center coordinates.

    throat_vector : string
        Dictionary key of the throat vector pointing along the length of the
        throats.

    throat_length : string
        Dictionary key of the throat length.

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    center = network[throat_centroid][throats]
    vector = network[throat_vector][throats]
    length = network[throat_length][throats]
    EP1 = center - 0.5 * length[:, _np.newaxis] * vector
    EP2 = center + 0.5 * length[:, _np.newaxis] * vector
    return {'head': EP1, 'tail': EP2}
