import scipy as _sp
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
    L = _ctc(target=target, pore_diameter=pore_diameter)
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    unit_vec = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
    EP1 = xyz[cn[:, 0]] + 0.5 * D1[:, _sp.newaxis] * unit_vec
    EP2 = xyz[cn[:, 1]] - 0.5 * D2[:, _sp.newaxis] * unit_vec
    # Handle overlapping pores
    overlap = L - 0.5 * (D1+D2) < 0
    mask = (D1 >= D2) & overlap
    EP2[mask] = EP1[mask]
    mask = (D1 < D2) & overlap
    EP1[mask] = EP2[mask]
    return {'head': EP1, 'tail': EP2}


def spherical_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter'):
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

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    This model should not be applied to true 2D networks. Use `circular_pores`
    model instead.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = _ctc(target=target, pore_diameter=pore_diameter)
    unit_vec = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    L1 = _sp.zeros_like(L)
    L2 = _sp.zeros_like(L)
    # Handle the case where Dt > Dp
    mask = Dt > D1
    L1[mask] = 0.5 * D1[mask]
    L1[~mask] = _sp.sqrt(D1[~mask]**2 - Dt[~mask]**2) / 2
    mask = Dt > D2
    L2[mask] = 0.5 * D2[mask]
    L2[~mask] = _sp.sqrt(D2[~mask]**2 - Dt[~mask]**2) / 2
    # Find throat endpoints
    EP1 = xyz[cn[:, 0]] + L1[:, None] * unit_vec
    EP2 = xyz[cn[:, 1]] - L2[:, None] * unit_vec
    # Handle throats w/ overlapping pores
    L1 = (4*L**2 + D1**2 - D2**2) / (8*L)
    # L2 = (4*L**2 + D2**2 - D1**2) / (8*L)
    h = (2*_sp.sqrt(D1**2/4 - L1**2)).real
    overlap = L - 0.5 * (D1+D2) < 0
    mask = overlap & (Dt < h)
    EP1[mask] = EP2[mask] = (xyz[cn[:, 0]] + L1[:, None] * unit_vec)[mask]
    return {'head': EP1, 'tail': EP2}


def circular_pores(target, pore_diameter='pore.diameter',
                   throat_diameter='throat.diameter'):
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

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    This model should only be applied to ture 2D networks.

    """
    return spherical_pores(target=target, pore_diameter=pore_diameter,
                           throat_diameter=throat_diameter)
