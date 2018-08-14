import scipy as _sp
from openpnm.utils import logging as _logging
_logger = _logging.getLogger(__name__)


def ctc(target, pore_diameter='pore.diameter'):
    r"""
    Calculate throat length assuming point-like pores, i.e. center-to-center
    distance between pores.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))
    return E


def straight(target, pore_diameter='pore.diameter', L_negative=1e-12):
    r"""
    Calculate throat length assuming no overlap between the throat and adjacent
    pores.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values

    L_negative : float
        The default throat length to use when negative lengths are found.  The
        default is 1 pm.  To accept negative throat lengths, set this value to
        ``None``.
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    E = ctc(target, pore_diameter=pore_diameter)
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    value = E-(D1+D2)/2.
    if _sp.any(value < 0) and L_negative is not None:
        _logger.warn('Negative throat lengths are calculated. Arbitrary ' +
                     'positive length assigned: ' + str(L_negative))
        Ts = _sp.where(value < 0)[0]
        value[Ts] = L_negative
    return value


def spherical_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter', L_negative=1e-12):
    r"""
    Calculate throat length, asusming pores are perfect spheres.

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

    L_negative : float
        The default throat length to use when negative lengths are found.  The
        default is 1 nm.  To accept negative throat lengths, set this value to
        ``None``.
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    Dt = network[throat_diameter][throats]
    if _sp.any([_sp.isnan(D1), _sp.isnan(D2), _sp.isnan(Dt)]):
        _logger.warn('Found spanner throats (spans between 2 geometries).' +
                     ' Either the other geometry is not defined yet, or it' +
                     ' does not have pore diameter values yet. Run' +
                     ' regenerate_models() on both geometries to fix.')
    L = ctc(target, pore_diameter=pore_diameter)
    L1 = _sp.zeros_like(L)
    L2 = _sp.zeros_like(L)
    # Handle the case where throat diam > pore diam
    mask = Dt > D1
    L1[mask] = D1[mask]/2
    L1[~mask] = _sp.sqrt(D1[~mask]**2 - Dt[~mask]**2) / 2
    mask = Dt > D2
    L2[mask] = D2[mask]/2
    L2[~mask] = _sp.sqrt(D2[~mask]**2 - Dt[~mask]**2) / 2
    # Handle throats w/ overlapping pores
    L1_temp = (4*L**2+D1**2-D2**2) / (8*L)
    L2_temp = (4*L**2+D2**2-D1**2) / (8*L)
    h_temp = (2*_sp.sqrt(D1**2/4 - L1_temp**2)).real
    mask_overlap = ((L - (D1+D2)/2) < 0) & (Dt < h_temp)
    L1[mask_overlap] = L1_temp[mask_overlap]
    L2[mask_overlap] = L2_temp[mask_overlap]
    # Calculate throat length
    Lt = L - (L1+L2)
    if _sp.any([L1 < 0, L2 < 0, Lt < 0]) and L_negative is not None:
        _logger.warn('Negative throat lengths are calculated. Arbitrary' +
                     ' positive length assigned: ' + str(L_negative))
    # Remove negative lengths
    Lt[Lt <= 0] = L_negative
    return Lt


def truncated_pyramid(target, pore_diameter='pore.diameter',
                      throat_diameter='throat.diameter', L_negative=1e-12):
    r"""
    Calculate throat length, assuming pores are spheres that are truncated
    to pyramids.

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

    L_negative : float
        The default throat length to use when negative lengths are found.  The
        default is 1 nm.  To accept negative throat lengths, set this value to
        ``None``.
    """
    return spherical_pores(target, pore_diameter=pore_diameter,
                           throat_diameter=throat_diameter,
                           L_negative=L_negative)


def circular_pores(target, pore_diameter='pore.diameter',
                   throat_diameter='throat.diameter', L_negative=1e-12):
    r"""
    Calculate throat length, assuming pores are perfect circles.


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

    L_negative : float
        The default throat length to use when negative lengths are found.  The
        default is 1 nm.  To accept negative throat lengths, set this value to
        ``None``.
    """
    return spherical_pores(target, pore_diameter=pore_diameter,
                           throat_diameter=throat_diameter,
                           L_negative=L_negative)


def boundary(target, pore_diameter='pore.diameter',
             throat_length='throat.length'):
    r"""
    Take the pore radius for the bulk pore, the throat length for the throat
    and use a minimum value for the boundary pores
    """
    net = target.project.network
    throats = net.map_throats(throats=target.Ts, origin=target)
    conns = net['throat.conns'][throats]
    tl = target[throat_length]
    p_lens = 0.999*net[pore_diameter][conns]/2
    return {'pore1': p_lens[:, 0], 'throat': tl, 'pore2': p_lens[:, 1]}
