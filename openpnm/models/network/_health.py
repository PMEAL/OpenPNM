"""
Health Checks
-------------

"""
import logging
import numpy as np
import pandas as pd
logger = logging.getLogger(__name__)

__all__ = [
    'bidirectional_throats',
    'cluster_number',
    'cluster_size',
    'count_coincident_pores',
    'duplicate_throats',
    'find_coincident_pores',
    'isolated_pores',
    'reversed_throats',
    'looped_throats',
    'headless_throats',
]


def cluster_number(network):
    r"""
    Assign a cluster number to each pore
    """
    from scipy.sparse import csgraph as csg
    am = network.create_adjacency_matrix(fmt='coo', triu=True)
    N, Cs = csg.connected_components(am, directed=False)
    return Cs


def cluster_size(network, cluster=None):
    r"""
    Find the size of the cluster to which each pore belongs

    Parameters
    ----------
    network : dict
        The Network
    cluster : str, optional
        Dict key pointing to the array containing the cluster number of each
        pore.  If not provided then it will be calculated.

    Returns
    -------
    cluster_size : ndarray
        An Np-long array containing the size of the cluster to which each pore
        belongs

    """
    if cluster is None:
        from scipy.sparse import csgraph as csg
        am = network.create_adjacency_matrix(fmt='coo', triu=True)
        N, cluster_num = csg.connected_components(am, directed=False)
    else:
        cluster_num = network[cluster]
    Cs, ind, N = np.unique(cluster_num, return_inverse=True, return_counts=True)
    values = N[ind]
    return values


def isolated_pores(network):
    r"""
    Find which pores, if any, are not connected to a throat
    """
    values = np.ones(network.Np, dtype=bool)
    hits = np.unique(network.conns)
    if np.any(hits >= network.Np):
        logger.warning("Some throats point to non-existent pores")
        hits = hits[hits < network.Np]
    values[hits] = False
    return values


def reversed_throats(network):
    r"""
    Find any throat connections that are pointing from j -> i where j > i
    """
    hits = network.conns[:, 0] > network.conns[:, 1]
    return hits


def looped_throats(network):
    r"""
    Find any throats that are connected to the same pore on both ends
    """
    hits = network.conns[:, 0] == network.conns[:, 1]
    return hits


def headless_throats(network):
    r"""
    Find any throats that point to a non-existent pore
    """
    hits = np.any(network.conns > (network.Np - 1), axis=1)
    return hits


def duplicate_throats(network):
    r"""
    Find repeat occurrences of throat connections
    """
    return pd.DataFrame(network.conns).duplicated().to_numpy()


def count_coincident_pores(network, thresh=1e-6):
    r"""
    Count number of pores that are spatially coincident with other pores

    Parameters
    ----------
    network : dict
        The Network
    thresh : float
        The distance below which two pores are considered spatially coincident

    Returns
    -------
    count : ndarray
        A numpy array of Np length containing the number of coincident pores
    """
    # This needs to be a bit complicated because it cannot be assumed
    # the coincident pores are topologically connected
    import scipy.spatial as sptl
    coords = network.coords
    tree = sptl.KDTree(coords)
    hits = tree.query_pairs(r=thresh)
    arr = np.array(list(hits)).flatten()
    v, n = np.unique(arr, return_counts=True)
    values = np.zeros(network.Np, dtype=int)
    values[v.astype(int)] = n
    return values


def find_coincident_pores(network, thresh=1e-6):
    r"""
    Find the indices of coincident pores

    Parameters
    ----------
    network : dict
        The Network
    thresh : float
        The distance below which two pores are considered spatially coincident

    Returns
    -------
    indices : list of lists
        One row corresponding to each pore, with each row listing the indices
        of any coincident pores.  An empty list means no pores were found
        within a distance of ``thresh``.

    """
    # This needs to be a bit complicated because it cannot be assumed
    # the coincident pores are topologically connected
    import scipy.spatial as sptl
    coords = network['pore.coords']
    tree = sptl.KDTree(coords)
    a = tree.sparse_distance_matrix(tree, max_distance=thresh,
                                    output_type='coo_matrix')
    a.data += 1.0
    a.setdiag(0)
    a.eliminate_zeros()
    a.data -= 1.0
    a = a.tolil()
    return a.rows


def bidirectional_throats(network):
    biTs = network['throat.conns'][:, 0] > network['throat.conns'][:, 1]
    return biTs
