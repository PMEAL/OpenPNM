import numpy as np
import scipy.stats as spst
import scipy.sparse as sprs
from scipy.sparse import csgraph
from collections import namedtuple


__all__ = [
    'trim_disconnected_clusters',
    'site_percolation_orig',
    'bond_percolation_orig',
    'ispercolating',
    'remove_isolated_clusters',
]


def trim_disconnected_clusters(b_labels, s_labels, inlets):
    r"""
    Computes actual node and edge occupancy based on connectivity to the given
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
    inlets : ndarray
        An array containing node indices that are to be treated as inlets.
        Any clusters labels not found in these nodes will be considered
        disconnected and set to -1.

    Returns
    -------
    occupancy : tuple of ndarrays
        The returned tuple contains boolean arrays of ``occupied_sites``
        and ``occupied_bonds``, after accounting for connection to the
        ``inlets``.

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


def remove_isolated_clusters(labels, inlets):
    r"""
    Finds cluster labels not attached to the inlets, and sets them to
    unoccupied (-1)

    Parameters
    ----------
    labels : tuple of node and edge labels
        This information is provided by the ``site_percolation`` or
        ``bond_percolation`` functions

    inlets : array_like
        A list of which nodes are inlets.  Can be a boolean mask or an
        array of indices.

    Returns
    -------
    A tuple containing a list of node and edge labels, with all clusters
    not connected to the inlets set to not occupied (-1).

    """
    # Identify clusters of invasion sites
    inv_clusters = np.unique(labels.sites[inlets])
    # Remove cluster numbers == -1, if any
    inv_clusters = inv_clusters[inv_clusters >= 0]
    # Find all pores in invading clusters
    p_invading = np.in1d(labels.sites, inv_clusters)
    labels.sites[~p_invading] = -1
    t_invading = np.in1d(labels.bonds, inv_clusters)
    labels.bonds[~t_invading] = -1
    return labels


def ispercolating(am, inlets, outlets, mode='site'):
    r"""
    Determines if a percolating clusters exists in the network spanning
    the given inlet and outlet nodes

    Parameters
    ----------
    am : adjacency_matrix
        The adjacency matrix with the ``data`` attribute indicating
        if an egde is occupied or not.
    inlets : array_like
        An array of indices indicating which nodes are part of the inlets
    outlets : array_like
        An array of indices indicating which nodes are part of the outlets
    mode : str
        Indicates which type of percolation to apply,. Options are:

        ===========  =====================================================
        mode         meaning
        ===========  =====================================================
        'site'       Applies site percolation
        'bond'       Applies bond percolation
        ===========  =====================================================

    """
    if am.format != 'coo':
        am = am.to_coo()
    ij = np.vstack((am.col, am.row)).T
    if mode.startswith('site'):
        occupied_sites = np.zeros(shape=am.shape[0], dtype=bool)
        occupied_sites[ij[am.data].flatten()] = True
        clusters = site_percolation_orig(ij, occupied_sites)
    elif mode.startswith('bond'):
        occupied_bonds = am.data
        clusters = bond_percolation_orig(ij, occupied_bonds)
    ins = np.unique(clusters.sites[inlets])
    if ins[0] == -1:
        ins = ins[1:]
    outs = np.unique(clusters.sites[outlets])
    if outs[0] == -1:
        outs = outs[1:]
    hits = np.in1d(ins, outs)
    return np.any(hits)


def bond_percolation_orig(conns, inds):
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
    inds : ndarray, boolean
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
    edge_mask = np.array(inds)
    Np = np.amax(conns) + 1
    # Find occupied sites based on occupied bonds
    # (the following 2 lines are not needed but worth keeping for future ref)
    # occupied_sites = np.zeros([Np, ], dtype=bool)
    # np.add.at(occupied_sites, ij[occupied_bonds].flatten(), True)
    adj_mat = sprs.csr_matrix((edge_mask, (conns[:, 0], conns[:, 1])),
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
    b_labels[~edge_mask] = -1
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


def site_percolation_orig(conns, inds):
    r"""
    Calculates the node and edge occupancy status for a site percolation
    process given a list of potentially occupied nodes.

    Parameters
    ----------
    conns : array_like
        An N x 2 array of [site_A, site_B] connections.  If two connected
        nodes are both occupied then they are part of the same cluster, as is
        the edge connecting them.
    inds : ndarray
        A boolean array with ``True`` values indicating occupied nodes.

    Returns
    -------
    A tuple containing a list of node and edge labels indicating which
    cluster each belongs to. A value of -1 indicates unoccupied.

    Notes
    -----
    The ``connected_components`` function of ``scipy.sparse.csgraph`` gives
    ALL nodes a cluster number whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    node_mask = np.array(inds)
    Np = np.size(inds)
    occupied_bonds = np.all(node_mask[conns], axis=1)
    adj_mat = sprs.csr_matrix((occupied_bonds, (conns[:, 0], conns[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    clusters[~node_mask] = -1
    s_labels = spst.rankdata(clusters + 1, method="dense") - 1
    if np.any(~node_mask):
        s_labels -= 1
    b_labels = np.amin(s_labels[conns], axis=1)
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)
