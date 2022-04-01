import scipy.sparse as sprs
import numpy as np
from scipy.sparse import csgraph


__all__ = [
    'bond_percolation',
    'find_connected_clusters',
    'find_trapped_clusters',
]


def bond_percolation(conns, occupied_bonds):
    r"""
    Calculates the site and bond occupancy status for a bond percolation
    process given a list of occupied bonds.

    Parameters
    ----------
    conns : array_like
        An N x 2 array of network connections, whether occupied or not.
        If the adjacency matrix is in the COO sparse format then this
        can be obtained using ``np.vstack(am.row, am.col).T``.  The
        matrix can be symmetric but only the upper triangular portion
        is kept.
    occupied_bonds: array containing bools
        A list of boolean values indicating whether a bond is occupied
        (``True``) or not (``False``). This must be the same length as
        ``conns``, but the lower triangular protion is ignored and the
        process is treated as undirected.

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to.  A value of -1 indicates uninvaded.

    Notes
    -----
    The ``connected_components`` function of ``scipy.sparse.csgraph`` will give
     a cluster number to ALL sites whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    from collections import namedtuple
    Np = np.amax(conns) + 1
    # Find occupied sites based on occupied bonds
    # (the following 2 lines are not needed but worth keeping for future ref)
    # occupied_sites = np.zeros([Np, ], dtype=bool)
    # np.add.at(occupied_sites, conns[occupied_bonds].flatten(), True)
    adj_mat = sprs.csr_matrix((occupied_bonds, (conns[:, 0], conns[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    # Clusters of size 1 only occur if all a site's bonds are uninvaded
    valid_clusters = np.bincount(clusters) > 1
    mapping = -np.ones(shape=(clusters.max()+1, ), dtype=int)
    mapping[valid_clusters] = np.arange(0, valid_clusters.sum())
    s_labels = mapping[clusters]
    # Bond inherit the cluster number of its connected sites
    b_labels = np.amin(s_labels[conns], axis=1)
    # Set bond cluster to -1 if not actually occupied
    b_labels[~occupied_bonds] = -1
    tup = namedtuple('cluster_labels', ('site_labels', 'bond_labels'))
    return tup(s_labels, b_labels)


def site_percolation(conns, occupied_sites):
    r"""
    Calculates the site and bond occupancy status for a site percolation
    process given a list of occupied sites

    Parameters
    ----------
    conns : array_like
        An N x 2 array of [site_A, site_B] connections.  If two connected
        sites are both occupied they are part of the same cluster, as it
        the bond connecting them.

    occupied_sites : bool
        A list indicating whether sites are occupied or not

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to.  A value of -1 indicates unoccupied.

    Notes
    -----
    The ``connected_components`` function of scipy.sparse.csgraph will give ALL
    sites a cluster number whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    from collections import namedtuple
    import scipy.stats as spst

    Np = np.size(occupied_sites)
    occupied_bonds = np.all(occupied_sites[conns], axis=1)
    adj_mat = sprs.csr_matrix((occupied_bonds, (conns[:, 0], conns[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    clusters[~occupied_sites] = -1
    s_labels = spst.rankdata(clusters + 1, method="dense") - 1
    if np.any(~occupied_sites):
        s_labels -= 1
    b_labels = np.amin(s_labels[conns], axis=1)
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


def find_connected_clusters(bond_labels, site_labels, inlets, asmask=True):
    hits = np.unique(site_labels[inlets])
    hits = hits[hits >= 0]
    occupied_bonds = np.isin(bond_labels, hits)
    occupied_sites = np.isin(site_labels, hits)
    if not asmask:
        occupied_bonds = occupied_bonds*(bond_labels + 1) - 1
        occupied_sites = occupied_sites*(site_labels + 1) - 1
    return occupied_sites, occupied_bonds


def find_trapped_clusters(conns, occupied_bonds, outlets):
    s_labels, b_labels = bond_percolation(conns, ~occupied_bonds)
    s_labels2, b_labels2 = find_connected_clusters(b_labels, s_labels,
                                                   outlets, asmask=False)
    # Set sites and bonds connected to outlets to -1, keeping
    s_labels[s_labels2 >= 0] = -1
    b_labels[b_labels2 >= 0] = -1
    return s_labels, b_labels













































