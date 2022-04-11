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


def find_path(network, pore_pairs, weights=None):
    am = network.create_adjacency_matrix(weights=weights)
    return queries.find_path(am, pairs=pore_pairs)


find_path.__doc__ = queries.find_path.__doc__


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


def trim_disconnected_clusters(**kwargs):
    return simulations.trim_disconnected_clusters(**kwargs)


trim_disconnected_clusters.__doc__ = \
    simulations.trim_disconnected_clusters.__doc__


def find_clusters(network, mask=[]):
    r"""
    Identify connected clusters of pores in the network.

    This method also returns a list of throat cluster numbers, which
    correspond to the cluster numbers of the pores to which the throat is
    connected.  Either site and bond percolation can be considered, see
    description of input arguments for details.

    Parameters
    ----------
    network : GenericNetwork
        The network
    mask : array_like, boolean
        A list of open bonds or sites (throats or pores).  If the mask is
        Np long, then the method will perform a site percolation to identify
        clusters, and if the mask is Nt long bond percolation will be
        performed.

    Returns
    -------
    p_labels, t_labels : tuple of ndarrays
        A tuple containing an Np long arra of pore cluster labels, and an
        Nt-long array of throat cluster labels. The label numbers correspond
        such that pores and throats with the same label are part of the same
        cluster.

    """
    # Parse the input arguments
    mask = np.array(mask, ndmin=1)
    if mask.dtype != bool:
        raise Exception('Mask must be a boolean array of Np or Nt length')
    # If pore mask was given perform site percolation
    if np.size(mask) == network.Np:
        (p_clusters, t_clusters) = \
            simulations.site_percolation_orig(network.conns, mask)
    # If pore mask was given perform bond percolation
    elif np.size(mask) == network.Nt:
        (p_clusters, t_clusters) = \
            simulations.bond_percolation_orig(network.conns, mask)
    else:
        raise Exception('Mask received was neither Nt nor Np long')
    return (p_clusters, t_clusters)
