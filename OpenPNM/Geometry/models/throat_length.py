r"""
===============================================================================
Submodule -- throat_length
===============================================================================

"""
import scipy as _sp
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


def straight(network, geometry, pore_diameter='pore.diameter',
             L_negative=1e-9, **kwargs):
    r"""
    Calculate throat length

    Parameters
    ----------
    L_negative : float
        The default throat length to use when negative lengths are found.  The
        default is 1 nm.  To accept negative throat lengths, set this value to
        ``None``.
    """
    # Initialize throat_property['length']
    throats = network.throats(geometry.name)
    pore1 = network['throat.conns'][:, 0]
    pore2 = network['throat.conns'][:, 1]
    C1 = network['pore.coords'][pore1]
    C2 = network['pore.coords'][pore2]
    E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))  # Euclidean distance between pores
    D1 = network[pore_diameter][pore1]
    D2 = network[pore_diameter][pore2]
    value = E-(D1+D2)/2.
    value = value[throats]
    if _sp.any(value < 0) and L_negative is not None:
        logger.warn('Negative throat lengths are calculated. Arbitrary positive' +
                    ' length assigned: ' + str(L_negative))
        Ts = _sp.where(value < 0)[0]
        value[Ts] = L_negative
    return value


def c2c(network, geometry, pore_centroid='pore.centroid',
        throat_centroid='throat.centroid', **kwargs):
    r"""
    Calculate the centre to centre distance from centroid of pore1 to centroid
    of throat to centroid of pore2.

    Parameters
    ----------
    pore_centroid and throat_centroid : string
        The dictionary keys to the arrays containing the pore and throat
        centroid coordinates.

    Notes
    -----
    This is tricky as connections are defined at Network level but centroids
    are stored on geometry. The pore and throat map relates the geometry index
    to the network index but we must look up the index of the map to go back
    to geometry index of the connected pores. This will probably break down
    when a throat connects two different geometries.
    """
    throats = geometry.map_throats(network, geometry.throats())
    connections = network['throat.conns'][throats]
    net_pore1 = connections[:, 0]
    net_pore2 = connections[:, 1]
    pore_centroids = network[pore_centroid]
    throat_centroids = network[throat_centroid][throats]
    v1 = throat_centroids-pore_centroids[net_pore1]
    v2 = throat_centroids-pore_centroids[net_pore2]
    value = _sp.ndarray(len(connections))
    for i in range(len(connections)):
        value[i] = _sp.linalg.norm(v1[i])+_sp.linalg.norm(v2[i])
    return value
