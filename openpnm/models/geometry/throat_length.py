r"""

.. autofunction:: openpnm.models.geometry.throat_length.piecewise
.. autofunction:: openpnm.models.geometry.throat_length.ctc
.. autofunction:: openpnm.models.geometry.throat_length.conduit_lengths

"""
import numpy as _np
from numpy.linalg import norm as _norm
from openpnm.utils import logging as _logging
_logger = _logging.getLogger(__name__)


def ctc(target):
    r"""
    Calculate throat length assuming point-like pores, i.e. center-to-center
    distance between pores. Also, this model assumes that pores and throat
    centroids are colinear.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    Returns
    -------
    value : NumPy ndarray
        Array containing throat length values.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    value = _norm(C1 - C2, axis=1)
    return value


def piecewise(target, throat_endpoints='throat.endpoints',
              throat_centroid='throat.centroid'):
    r"""
    Calculate throat length from end points and optionally a centroid

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_endpoints : string
        Dictionary key of the throat endpoint values.

    throat_centroid : string
        Dictionary key of the throat centroid values, optional.

    Returns
    -------
    Lt : ndarray
        Array containing throat lengths for the given geometry.

    Notes
    -----
    (1) By default, the model assumes that the centroids of pores and the
    connecting throat in each conduit are colinear.

    (2) If `throat_centroid` is passed, the model accounts for the extra
    length. This could be useful for Voronoi or extracted networks.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    # Get throat endpoints
    EP1 = network[throat_endpoints + '.head'][throats]
    EP2 = network[throat_endpoints + '.tail'][throats]
    # Calculate throat length
    Lt = _norm(EP1 - EP2, axis=1)
    # Handle the case where pores & throat centroids are not colinear
    try:
        Ct = network[throat_centroid][throats]
        Lt = _norm(Ct - EP1, axis=1) + _norm(Ct - EP2, axis=1)
    except KeyError:
        pass
    return Lt


def conduit_lengths(target, throat_endpoints='throat.endpoints',
                    throat_length='throat.length',
                    throat_centroid='throat.centroid'):
    r"""
    Calculate conduit lengths. A conduit is defined as half pore + throat
    + half pore.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_endpoints : string
        Dictionary key of the throat endpoint values.

    throat_diameter : string
        Dictionary key of the throat length values.

    throat_length : string (optional)
        Dictionary key of the throat length values.  If not given then the
        direct distance bewteen the two throat end points is used.

    Returns
    -------
    Dictionary containing conduit lengths, which can be accessed via the dict
    keys 'pore1', 'pore2', and 'throat'.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    # Get pore coordinates
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    # Get throat endpoints and length
    EP1 = network[throat_endpoints + '.head'][throats]
    EP2 = network[throat_endpoints + '.tail'][throats]
    try:
        # Look up throat length if given
        Lt = network[throat_length][throats]
    except KeyError:
        # Calculate throat length otherwise based on piecewise model
        Lt = piecewise(target, throat_endpoints, throat_centroid)
    # Calculate conduit lengths for pore 1 and pore 2
    L1 = _norm(C1 - EP1, axis=1)
    L2 = _norm(C2 - EP2, axis=1)
    return {'pore1': L1, 'throat': Lt, 'pore2': L2}


def classic(target, pore_diameter='pore.diameter'):
    r"""
    Find throat length as the pore-to-pore center distance, less the radii of
    each pore.

    Parameters
    ----------
    target : OpenPNM object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    ctc_dist = ctc(target)
    value = ctc_dist - network[pore_diameter][cn].sum(axis=1)/2
    return value
