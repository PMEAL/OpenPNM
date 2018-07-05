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
    throats = target.map_throats(target['throat._id'])
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))
    return E


def straight(target, pore_diameter='pore.diameter', L_negative=1e-12):
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
    throats = target.map_throats(target['throat._id'])
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
                    throat_diameter='throat.diameter'):
    r"""
    Calculate conduit lengths, i.e. pore 1 length, throat length,
    and pore 2 length, assuming that pores are spheres.

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

    """
    network = target.project.network
    throats = target.map_throats(target['throat._id'])
    cn = network['throat.conns'][throats]
    d1 = network[pore_diameter][cn[:, 0]]
    d2 = network[pore_diameter][cn[:, 1]]
    dt = network[throat_diameter][throats]
    L = ctc(target, pore_diameter=pore_diameter)
    L1 = _sp.sqrt(d1**2 - dt**2) / 2            # Effective length of pore 1
    L2 = _sp.sqrt(d2**2 - dt**2) / 2            # Effective length of pore 2
    Lt = L - (L1+L2)                            # Effective length of throat
    # Handling throats w/ overlapping pores
    L1temp = (4*L**2+d1**2-d2**2) / (8*L)
    L2temp = (4*L**2+d2**2-d1**2) / (8*L)
    htemp = (2*_sp.sqrt(d1**2/4 - L1temp**2)).real
    mask_overlap = ((L - (d1+d2)/2) < 0) & (dt < htemp)
    L1[mask_overlap] = L1temp[mask_overlap]
    L2[mask_overlap] = L2temp[mask_overlap]
    Lt[mask_overlap] = 1e-12
    return {'pore1': L1, 'throat': Lt, 'pore2': L2}


def truncated_pyramid(target, pore_diameter='pore.diameter',throat_diameter='throat.diameter'):
    r"""
    Calculate conduit lengths, i.e. pore 1 length, throat length,
    and pore 2 length, assuming that pores are pyramid.

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

    """
    return spherical_pores(target, pore_diameter=pore_diameter,
                           throat_diameter=throat_diameter)


def circular_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter'):
    r"""
    Calculate conduit lengths, i.e. pore 1 length, throat length,
    and pore 2 length, assuming that pores are circles.

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

    """
    return spherical_pores(target, pore_diameter=pore_diameter,
                           throat_diameter=throat_diameter)
