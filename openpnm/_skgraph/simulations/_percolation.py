from collections import namedtuple

import numpy as np
import scipy.sparse as sprs
from scipy.sparse import csgraph

from openpnm._skgraph.operations import split_edges

__all__ = [
    'trim_disconnected_clusters',
    'ispercolating',
    'remove_isolated_clusters',
    'bond_percolation',
    'site_percolation',
    # 'find_connected_clusters',
    # 'find_trapped_bonds',
    # 'find_trapped_sites',
    'find_connected_clusters',
]


def bond_percolation(conns, occupied_bonds):
    r"""
    Assigns cluster numbers to sites and bonds acccording to a bond
    percolation process, given a list of occupied bonds.

    Parameters
    ----------
    conns : array_like
        An N x 2 array connections. Any sites connected to an occupied bond
        will also be considered occupied and given the same cluster number
        as the bond.
    occupied_bonds : ndarray
        A boolean array with one element for each bond, with ``True`` values
        indicating that a bond is occupied

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to. A value of -1 indicates uninvaded.

    Notes
    -----
    The ``connected_components`` function of ``scipy.sparse.csgraph`` will give
    a cluster number to ALL bonds whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    Np = np.amax(conns) + 1
    # Find occupied sites based on status of shared bonds
    occupied_sites = np.zeros([Np, ], dtype=bool)
    occupied_sites[conns[occupied_bonds].flatten()] = True
    # Perform cluster labeling of network
    adj_mat = sprs.csr_matrix((occupied_bonds.astype(int),
                               (conns[:, 0], conns[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    # Set cluster number of unoccupied sites to -1
    s_labels = (clusters + 1)*occupied_sites - 1
    # Bonds inherit the cluster number of its connected sites
    b_labels = np.amin(s_labels[conns], axis=1)
    # Set cluster number of unoccupied bonds to -1
    b_labels[~occupied_bonds] = -1
    tup = namedtuple('cluster_labels', ('site_labels', 'bond_labels'))
    return tup(s_labels, b_labels)


def site_percolation(conns, occupied_sites):
    r"""
    Assigns cluster numbers to sites and bonds acccording to a site
    percolation process, given a list of occupied sites.

    Parameters
    ----------
    conns : array_like
        An N x 2 array connections. If two connected sites are both occupied
        they are part of the same cluster, as is the bond connecting them.
    occupied_sites : ndarray
        A boolean array with one element for each site, with ``True`` values
        indicating that a site is occupied

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to.  A value of -1 indicates unoccupied.

    Notes
    -----
    The ``connected_components`` function of ``scipy.sparse.csgraph`` will
    give ALL sites a cluster number whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    Np = np.size(occupied_sites)
    # Find bond occupancy based on status of its connected sites
    occupied_bonds = np.all(occupied_sites[conns], axis=1)
    # Perform cluster labeling of network
    adj_mat = sprs.csr_matrix((occupied_bonds, (conns[:, 0], conns[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    # Set cluster number of unoccupied sites to -1
    s_labels = (clusters + 1)*occupied_sites - 1
    # Bonds inherit the cluster number of its connected sites
    b_labels = np.amin(s_labels[conns], axis=1)
    # Set cluster number of unoccupied bonds to -1
    b_labels[~occupied_bonds] = -1
    tup = namedtuple('cluster_labels', ('site_labels', 'bond_labels'))
    return tup(s_labels, b_labels)


def mixed_percolation(
    conns,
    occupied_sites,
    occupied_bonds
):  # pragma: no cover
    new_conns = split_edges(conns)[0]
    new_sites = np.hstack((occupied_sites, occupied_bonds))
    s, b = site_percolation(conns=new_conns, occupied_sites=new_sites)
    s_labels = s[:occupied_sites.shape[0]]
    b_labels = s[occupied_sites.shape[0]:]
    return s_labels, b_labels


def find_connected_clusters(
    bond_labels,
    site_labels,
    inlets,
    asmask=True
):  # pragma: no cover
    hits = np.unique(site_labels[inlets])
    hits = hits[hits >= 0]
    occupied_bonds = np.isin(bond_labels, hits)
    occupied_sites = np.isin(site_labels, hits)
    if not asmask:
        occupied_bonds = occupied_bonds*(bond_labels + 1) - 1
        occupied_sites = occupied_sites*(site_labels + 1) - 1
    return occupied_sites, occupied_bonds


def find_trapped_bonds(conns, outlets, occupied_bonds):  # pragma: no cover
    s_labels, b_labels = bond_percolation(conns, ~occupied_bonds)
    s_labels2, b_labels2 = find_connected_clusters(b_labels, s_labels,
                                                   outlets, asmask=False)
    # Set sites and bonds connected to outlets to -1, keeping
    s_labels[s_labels2 >= 0] = -1
    b_labels[b_labels2 >= 0] = -1
    return s_labels, b_labels


def find_trapped_sites(conns, outlets, occupied_sites):  # pragma: no cover
    s_labels, b_labels = site_percolation(conns, ~occupied_sites)
    s_labels2, b_labels2 = find_connected_clusters(b_labels, s_labels,
                                                   outlets, asmask=False)
    # Set sites and bonds connected to outlets to -1, keeping
    s_labels[s_labels2 >= 0] = -1
    b_labels[b_labels2 >= 0] = -1
    return s_labels, b_labels


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
        The returned tuple containing arrays of cluster numbers of
        ``occupied_sites`` and ``occupied_bonds``, after accounting for
        connection to the ``inlets``.

    Notes
    -----
    The ``b_labels`` and ``s_labels`` arrays are returned from the
    ``bond_percolation`` or ``site_percolation`` function.

    """
    hits = np.unique(s_labels[inlets])
    hits = hits[hits >= 0]
    occupied_bonds = np.isin(b_labels, hits)*(b_labels + 1) - 1
    occupied_sites = np.isin(s_labels, hits)*(s_labels + 1) - 1
    r = namedtuple('cluster_labels', ('site_labels', 'bond_labels'))
    return r(occupied_sites, occupied_bonds)


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
    inv_clusters = np.unique(labels.site_labels[inlets])
    # Remove cluster numbers == -1, if any
    inv_clusters = inv_clusters[inv_clusters >= 0]
    # Find all pores in invading clusters
    p_invading = np.in1d(labels.site_labels, inv_clusters)
    labels.site_labels[~p_invading] = -1
    t_invading = np.in1d(labels.bond_labels, inv_clusters)
    labels.bond_labels[~t_invading] = -1
    return labels


def ispercolating(conns, occupied, inlets, outlets):
    r"""
    Determines if a percolating cluster exists in the network spanning
    the given inlet and outlet nodes

    Parameters
    ----------
    conns : array_like
        An N x 2 array connections. If two connected sites are both occupied
        they are part of the same cluster, as is the bond connecting them.
    occupied : array_like
        A boolean array with ``True`` values indicating if a bond or site
        is occupied.  If the length of this array is equal to the number
        of bonds (i.e. ``conns.shape[0]``) then bond percolation is assumed,
        otherwise site percolation is assumed.
    inlets : array_like
        An array of indices indicating which nodes are part of the inlets
    outlets : array_like
        An array of indices indicating which nodes are part of the outlets
    """
    if occupied.size == conns.shape[0]:
        mode = 'bond'
    else:
        mode = 'site'
    if mode.startswith('site'):
        clusters = site_percolation(conns, occupied)
    elif mode.startswith('bond'):
        clusters = bond_percolation(conns, occupied)
    ins = np.unique(clusters.site_labels[inlets])
    if ins[0] == -1:
        ins = ins[1:]
    outs = np.unique(clusters.site_labels[outlets])
    if outs[0] == -1:
        outs = outs[1:]
    hits = np.in1d(ins, outs)
    return np.any(hits)
