import scipy.sparse as sprs
import numpy as np
from scipy.sparse import csgraph


__all__ = [
    'bond_percolation',
    'trim_disconnected_clusters',
]


def bond_percolation(conns, occupied_bonds):
    r"""
    Calculates the site and bond occupancy status for a bond percolation
    process given a list of occupied bonds.

    Parameters
    ----------
    conns : array_like
        An N x 2 array of connections whether occupied or not.
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
    # np.add.at(occupied_sites, ij[occupied_bonds].flatten(), True)
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
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


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
        unoccupied.

    Returns
    -------
    occupancy : tuple of ndarrays
        The return tuple contains boolean arrays of ``occupied_sites``
        and ``occupied_bonds``, after accounting for connection to the
        inlets.

    Notes
    -----
    The ``b_labels`` and ``s_labels`` arrays are returned from the
    ``bond_percolation`` and ``site_percolation`` functions.

    """
    hits = np.unique(s_labels[inlets])
    hits = hits[hits >= 0]
    occupied_bonds = np.isin(b_labels, hits)
    occupied_sites = np.isin(s_labels, hits)
    return occupied_sites, occupied_bonds






















