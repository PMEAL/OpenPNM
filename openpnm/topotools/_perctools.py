import logging
import numpy as np
import scipy.sparse as sprs
from scipy.sparse import csgraph
from openpnm.utils import PrintableDict, Workspace
from openpnm._skgraph import simulations
from openpnm._skgraph import queries
from collections import namedtuple


logger = logging.getLogger(__name__)
ws = Workspace()
__all__ = [
    'ispercolating',
    'find_isolated_clusters',
    'trim_disconnected_clusters',
    'site_percolation',
    'bond_percolation',
    'find_clusters',
    'find_path',
]


def find_path(network, pore_pairs, weights=None):
    return queries.find_path(network=network, pairs=pore_pairs, weights=weights)


find_path.__doc__ = queries.find_path.__doc__


def ispercolating(network, inlets, outlets):
    if np.array(inlets).dtype == bool:
        inlets = np.where(inlets)[0]
    if np.array(outlets).dtype == bool:
        outlets = np.where(outlets)[0]
    flag = simulations.ispercolating(
        conns=network.conns,
        occupied=np.ones(network.Nt, dtype=bool),
        inlets=inlets,
        outlets=outlets
    )
    return flag


ispercolating.__doc__ = simulations.ispercolating.__doc__


def site_percolation(network, occupied_sites):
    return simulations.site_percolation(network.conns, occupied_sites)


site_percolation.__doc__ = simulations.site_percolation.__doc__


def bond_percolation(network, occupied_bonds):
    return simulations.bond_percolation(network.conns, occupied_bonds)


bond_percolation.__doc__ = simulations.bond_percolation.__doc__


def trim_disconnected_clusters(**kwargs):
    return simulations.trim_disconnected_clusters(**kwargs)


trim_disconnected_clusters.__doc__ = \
    simulations.trim_disconnected_clusters.__doc__


def find_isolated_clusters(network, mask, inlets):
    r"""
    Identifies pores and throats that are invaded but not connected to the inlets

    Parameters
    ----------
    network : dict
        The OpenPNM Network
    mask : ndarray
        A boolean mask of either Nt or Np length with ``True`` values
        indicating invaded bonds or sites. If this array is Nt-long then
        then bond percolation is used to identify clusters, whereas site
        percolation is used if it is Np-long.
    inlets : ndarray
        A array containing indices of the pores which define the inlets. Any
        clusters not connected to these sites are considered isolated.

    Returns
    -------
    sites : ndarray
        An ndarray containing the indices of invaded pores which are not
        connected to the given ``inlets``.
    """
    labels = find_clusters(network=network, mask=mask)
    isolated = np.in1d(labels.pore_labels, labels.pore_labels[inlets], invert=True)
    isolated = np.where(isolated)[0]
    return isolated


def find_clusters(network, mask=[]):
    r"""
    Identify connected clusters of pores and throats in the network.

    Either site and bond percolation can be considered, see description of
    ``mask`` argument for details.

    Parameters
    ----------
    network : Network
        The network
    mask : array_like, boolean
        A list of open bonds or sites (throats or pores). If the mask is
        Np-long, then the method will perform a site percolation to identify
        clusters, and if the mask is Nt-long bond percolation will be
        performed.

    Returns
    -------
    p_labels, t_labels : tuple of ndarrays
        A tuple containing an Np-long array of pore cluster labels, and an
        Nt-long array of throat cluster labels. The label numbers correspond
        such that pores and throats with the same label are part of the same
        cluster. Uninvaded locations are set to -1.

    """
    # Parse the input arguments
    mask = np.array(mask, ndmin=1)
    if mask.dtype != bool:
        raise Exception('Mask must be a boolean array of Np or Nt length')
    # If pore mask was given perform site percolation
    if np.size(mask) == network.Np:
        (p_clusters, t_clusters) = \
            simulations.site_percolation(network.conns, mask)
    # If pore mask was given perform bond percolation
    elif np.size(mask) == network.Nt:
        (p_clusters, t_clusters) = \
            simulations.bond_percolation(network.conns, mask)
    else:
        raise Exception('Mask received was neither Nt nor Np long')
    result = namedtuple('result', ('pore_labels', 'throat_labels'))
    result.pore_labels = p_clusters
    result.throat_labels = t_clusters
    return result
