import logging
import numpy as np
import scipy.sparse as sprs
from scipy.sparse import csgraph
from openpnm.utils import PrintableDict, Workspace
from openpnm._skgraph import simulations
from openpnm._skgraph import queries


logger = logging.getLogger(__name__)
ws = Workspace()
__all__ = [
    'ispercolating',
    'remove_isolated_clusters',
    'site_percolation',
    'bond_percolation',
    'find_clusters',
    'find_path',
]


def ispercolating(**kwargs):
    return simulations.ispercolating(**kwargs)


ispercolating.__doc__ = simulations.ispercolating.__doc__


def remove_isolated_clusters(**kwargs):
    return simulations.remove_isolated_clusters(**kwargs)


remove_isolated_clusters.__doc__ = simulations.remove_isolated_clusters.__doc__


def site_percolation(ij, occupied_sites):
    return simulations.site_percolation_orig(conns=ij, inds=occupied_sites)


site_percolation.__doc__ = simulations.site_percolation_orig.__doc__


def bond_percolation(conns, occupied_bonds):
    return simulations.bond_percolation_orig(conns=conns, inds=occupied_bonds)


bond_percolation.__doc__ = simulations.bond_percolation_orig.__doc__


def trim_disconnected_clusters(b_labels, s_labels, inlets):
    r"""
    Computes actual site and bond occupancy based on connectivity to the given
    inlets

    Parameters
    ----------
    b_labels : ndarray
        An array of cluster labels assigned to each bond.  -1 indicates
        unoccupied
    s_labels : ndarray
        An array of cluster labels assigned to each site. -1 indicates
        unoccupied. Site cluster numbers must correspond to the bond
        clusters, such that if bond j has a cluster number N, then both
        sites on each end of j are also labeled N.

    Returns
    -------
    occupancy : tuple of ndarrays
        The returned tuple contains boolean arrays of ``occupied_sites``
        and ``occupied_bonds``, after accounting for connection to the
        inlets.

    Notes
    -----
    The ``b_labels`` and ``s_labels`` arrays are returned from the
    ``bond_percolation`` function.

    """
    hits = np.unique(s_labels[inlets])
    hits = hits[hits >= 0]
    occupied_bonds = np.isin(b_labels, hits)
    occupied_sites = np.isin(s_labels, hits)
    return occupied_sites, occupied_bonds


def find_clusters(network, mask=[], t_labels=False):
    r"""
    Identify connected clusters of pores in the network.  This method can
    also return a list of throat cluster numbers, which correspond to the
    cluster numbers of the pores to which the throat is connected.  Either
    site and bond percolation can be considered, see description of input
    arguments for details.

    Parameters
    ----------
    network : GenericNetwork
        The network

    mask : array_like, boolean
        A list of active bonds or sites (throats or pores).  If the mask is
        Np long, then the method will perform a site percolation, and if
        the mask is Nt long bond percolation will be performed.

    Returns
    -------
    A tuple containing an Np long list of pore cluster labels, and an Nt-long
    list of throat cluster labels.  The label numbers correspond such that
    pores and throats with the same label are part of the same cluster.

    """
    # Parse the input arguments
    mask = np.array(mask, ndmin=1)
    if mask.dtype != bool:
        raise Exception('Mask must be a boolean array of Np or Nt length')

    # If pore mask was given perform site percolation
    if np.size(mask) == network.Np:
        (p_clusters, t_clusters) = _site_percolation(network, mask)
    # If pore mask was given perform bond percolation
    elif np.size(mask) == network.Nt:
        (p_clusters, t_clusters) = _bond_percolation(network, mask)
    else:
        raise Exception('Mask received was neither Nt nor Np long')

    return (p_clusters, t_clusters)


def _site_percolation(network, pmask):
    r"""
    This private method is called by 'find_clusters'
    """
    # Find throats that produce site percolation
    conns = np.copy(network['throat.conns'])
    conns[:, 0] = pmask[conns[:, 0]]
    conns[:, 1] = pmask[conns[:, 1]]
    # Only if both pores are True is the throat set to True
    tmask = np.all(conns, axis=1)

    # Perform the clustering using scipy.sparse.csgraph
    csr = network.create_adjacency_matrix(weights=tmask, fmt='csr',
                                          drop_zeros=True)
    clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                 directed=False)[1]

    # Adjust cluster numbers such that non-invaded pores are labelled -1
    # Note: The following line also takes care of assigning cluster numbers
    # to single isolated invaded pores
    p_clusters = (clusters + 1)*(pmask) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']]
    ind = (t_clusters[:, 0] == t_clusters[:, 1])
    t_clusters = t_clusters[:, 0]
    # Label non-invaded throats with -1
    t_clusters[~ind] = -1

    return (p_clusters, t_clusters)


def _bond_percolation(network, tmask):
    r"""
    This private method is called by 'find_clusters'
    """
    # Perform the clustering using scipy.sparse.csgraph
    csr = network.create_adjacency_matrix(weights=tmask, fmt='csr',
                                          drop_zeros=True)
    clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                 directed=False)[1]

    # Convert clusters to a more usable output:
    # Find pores attached to each invaded throats
    Ps = network.find_connected_pores(throats=tmask, flatten=True)
    # Adjust cluster numbers such that non-invaded pores are labelled -1
    p_clusters = (clusters + 1)*(network.to_mask(pores=Ps).astype(int)) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']][:, 0]
    # Label non-invaded throats with -1
    t_clusters[~tmask] = -1

    return (p_clusters, t_clusters)


def find_path(network, pore_pairs, weights=None):
    am = network.create_adjacency_matrix(weights=weights)
    return queries.find_path(am, pairs=pore_pairs)


find_path.__doc__ = queries.find_path.__doc__
