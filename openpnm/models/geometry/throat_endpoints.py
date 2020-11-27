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


def spherical_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter',
                    throat_centroid='throat.centroid'):
    r"""
    Calculate the coordinates of throat endpoints, assuming spherical pores.
    This model accounts for the overlapping lens between pores and throats.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.

    throat_centroid : string, optional
        Dictionary key of the throat centroid values. See the notes.

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    (1) This model should not be applied to true 2D networks. Use
    ``circular_pores`` model instead.

    (2) By default, this model assumes that throat centroid and pore
    coordinates are colinear. If that's not the case, such as in extracted
    networks, ``throat_centroid`` could be passed as an optional argument, and
    the model takes care of the rest.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = _ctc(target=target) + 1e-15
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    L1 = _np.zeros_like(L)
    L2 = _np.zeros_like(L)
    # Handle the case where Dt > Dp
    mask = Dt > D1
    L1[mask] = 0.5 * D1[mask]
    L1[~mask] = _np.sqrt(D1[~mask]**2 - Dt[~mask]**2) / 2
    mask = Dt > D2
    L2[mask] = 0.5 * D2[mask]
    L2[~mask] = _np.sqrt(D2[~mask]**2 - Dt[~mask]**2) / 2
    # Handle non-colinear pores and throat centroids
    try:
        TC = network[throat_centroid][throats]
        LP1T = _np.linalg.norm(TC - xyz[cn[:, 0]], axis=1) + 1e-15
        LP2T = _np.linalg.norm(TC - xyz[cn[:, 1]], axis=1) + 1e-15
        unit_vec_P1T = (TC - xyz[cn[:, 0]]) / LP1T[:, None]
        unit_vec_P2T = (TC - xyz[cn[:, 1]]) / LP2T[:, None]
    except KeyError:
        unit_vec_P1T = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
        unit_vec_P2T = -1 * unit_vec_P1T
    # Find throat endpoints
    EP1 = xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T
    EP2 = xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T
    # Handle throats w/ overlapping pores
    L1 = (4 * L**2 + D1**2 - D2**2) / (8 * L)
    L2 = (4 * L**2 + D2**2 - D1**2) / (8 * L)
    h = (2 * _np.sqrt(D1**2 / 4 - L1**2)).real
    overlap = L - 0.5 * (D1 + D2) < 0
    mask = overlap & (Dt < h)
    EP1[mask] = (xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T)[mask]
    EP2[mask] = (xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T)[mask]
    return {'head': EP1, 'tail': EP2}


def circular_pores(target, pore_diameter='pore.diameter',
                   throat_diameter='throat.diameter',
                   throat_centroid='throat.centroid'):
    r"""
    Calculate the coordinates of throat endpoints, assuming circular pores.
    This model accounts for the overlapping lens between pores and throats.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.

    throat_centroid : string, optional
        Dictionary key of the throat centroid values. See the notes.

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    (1) This model should only be applied to ture 2D networks.

    (2) By default, this model assumes that throat centroid and pore
    coordinates are colinear. If that's not the case, such as in extracted
    networks, `throat_centroid` could be passed as an optional argument, and
    the model takes care of the rest.

    """
    return spherical_pores(target, pore_diameter=pore_diameter,
                           throat_diameter=throat_diameter)


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
