import scipy as _sp
from openpnm.core import logging as _logging
_logger = _logging.getLogger(__name__)


def straight(target, pore_diameter='pore.diameter', L_negative=1e-12):
    r"""
    Calculate throat length

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    L_negative : float
        The default throat length to use when negative lengths are found.  The
        default is 1 nm.  To accept negative throat lengths, set this value to
        ``None``.
    """
    network = target.project.network
    # Initialize throat_property['length']
    throats = network.throats(target.name)
    pore1 = network['throat.conns'][:, 0]
    pore2 = network['throat.conns'][:, 1]
    C1 = network['pore.coords'][pore1]
    C2 = network['pore.coords'][pore2]
    E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))  # Euclidean distance
    D1 = network[pore_diameter][pore1]
    D2 = network[pore_diameter][pore2]
    value = E-(D1+D2)/2.
    value = value[throats]
    if _sp.any(value < 0) and L_negative is not None:
        _logger.warn('Negative throat lengths are calculated. Arbitrary ' +
                     'positive length assigned: ' + str(L_negative))
        Ts = _sp.where(value < 0)[0]
        value[Ts] = L_negative
    return value


def endpoints(target,
              throat_vector='throat.vector',
              pore_diameter='pore.diameter'):
    r"""
    Find the start and end point of each throat by offset each pore coordinate
    by it's radius along the throat vector.

    Parameters
    ----------
    throat_vector : string
        The normalized vector desribing the throat orientation.  The default
        is 'throat.vector', but if not found it will be calculated.

    pore_diameter : string
        The diameter of the pores.  The default is 'pore.diameter'.

    Examples
    --------
    The following code snippet can be used to visualize the throats as line
    segments:

    """
    net = target.project.network
    P12 = net['throat.conns']
    try:
        unit_vec = net[throat_vector]
    except KeyError:
        vec = _sp.diff(net['pore.coords'][P12], axis=1).squeeze()
        mag = _sp.atleast_2d(_sp.sqrt(_sp.sum(_sp.square(vec), axis=1))).T
        unit_vec = vec/mag
    # Find the principle angles of each vector
    rads = _sp.arccos(unit_vec)
    # Find endpoints of each throat by offsetting pore centers by their radii
    pdia1 = _sp.atleast_2d(net[pore_diameter][P12[:, 0]]).T
    pdia2 = _sp.atleast_2d(net[pore_diameter][P12[:, 1]]).T
    tends1 = net['pore.coords'][P12[:, 0]] + 0.5*pdia1*_sp.cos(rads)
    tends2 = net['pore.coords'][P12[:, 1]] - 0.5*pdia2*_sp.cos(rads)
    tends1 = _sp.swapaxes(_sp.atleast_3d(tends1), 2, 1)
    tends2 = _sp.swapaxes(_sp.atleast_3d(tends2), 2, 1)
    vals = _sp.concatenate((tends1, tends2), axis=1)
    return vals[net.throats(target.name)]


def conduit_lengths(target,
                    throat_endpoints='throat.endpoints',
                    throat_centroid='throat.centroid'):
    r"""
    Return the respective lengths of the conduit components in the form of a
    Nt x 3 array of [P1_len, T_len, P2_len]

    Parameters
    ----------
    throat_endpoints : string
        An Nt x 2 x 3 array containing [X1, Y1, Z1], [X2, Y2, Z2] pairs on each
        row indicating the start and end points of a throat.  These should be
        arranged in the same order as 'throat.conns', so if pore 1 is
        connected to pore 2, then this array should have the throat endpoint
        associated with pore 1 first, and with pore 2 second.

    throat_centroids : string
        An Nt x 3 array containing the [X, Y, Z] coordinates indicatding the
        center of the throat.  This value is only used if ``throat_endpoints``
        is not available, and will result a throat length of zero.

    """
    network = target.project.network
    P12 = network['throat.conns']
    pcoords = network['pore.coords']
    # Use endpoints to find pore lengths
    try:
        tends = target[throat_endpoints]
    except KeyError:
        tends = target[throat_centroid]
        tends = _sp.swapaxes(_sp.atleast_3d(tends), 2, 1)
        tends = _sp.concatenate((tends, tends), axis=1)
    plen1 = _sp.sqrt(_sp.sum(_sp.square(pcoords[P12[:, 0]]-tends[:, 0]), 1))
    plen2 = _sp.sqrt(_sp.sum(_sp.square(pcoords[P12[:, 1]]-tends[:, 1]), 1))
    # Use endpoints to find throat length
    tlen = _sp.sqrt(_sp.sum(_sp.square(_sp.diff(tends, axis=1)),
                            axis=2)).squeeze()
    return _sp.vstack((plen1, tlen, plen2)).T
