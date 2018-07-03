import scipy as _sp
from openpnm.core import logging as _logging
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
    throats = network.throats(target.name)
    cn = network['throat.conns']
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))
    return E[throats]


def straight(target, pore_diameter='pore.diameter', L_negative=1e-9):
    r"""
    Calculate throat length

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
