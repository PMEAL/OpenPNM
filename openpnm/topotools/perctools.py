import numpy as np
import scipy.sparse as sprs
from scipy.sparse import csgraph
from openpnm.utils import PrintableDict, logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


def ispercolating(am, inlets, outlets, mode='site'):
    r"""
    Determines if a percolating clusters exists in the network spanning
    the given inlet and outlet sites

    Parameters
    ----------
    am : adjacency_matrix
        The adjacency matrix with the ``data`` attribute indicating
        if a bond is occupied or not

    inlets : array_like
        An array of indices indicating which sites are part of the inlets

    outlets : array_like
        An array of indices indicating which sites are part of the outlets

    mode : string
        Indicates which type of percolation to apply, either `'site'` or
        `'bond'`

    """
    if am.format != 'coo':
        am = am.to_coo()
    ij = np.vstack((am.col, am.row)).T
    if mode.startswith('site'):
        occupied_sites = np.zeros(shape=am.shape[0], dtype=bool)
        occupied_sites[ij[am.data].flatten()] = True
        clusters = site_percolation(ij, occupied_sites)
    elif mode.startswith('bond'):
        occupied_bonds = am.data
        clusters = bond_percolation(ij, occupied_bonds)
    ins = np.unique(clusters.sites[inlets])
    if ins[0] == -1:
        ins = ins[1:]
    outs = np.unique(clusters.sites[outlets])
    if outs[0] == -1:
        outs = outs[1:]
    hits = np.in1d(ins, outs)
    return np.any(hits)


def remove_isolated_clusters(labels, inlets):
    r"""
    Finds cluster labels not attached to the inlets, and sets them to
    unoccupied (-1)

    Parameters
    ----------
    labels : tuple of site and bond labels
        This information is provided by the ``site_percolation`` or
        ``bond_percolation`` functions

    inlets : array_like
        A list of which sites are inlets.  Can be a boolean mask or an
        array of indices.

    Returns
    -------
    A tuple containing a list of site and bond labels, with all clusters
    not connected to the inlet sites set to not occupied.

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


def site_percolation(ij, occupied_sites):
    r"""
    Calculates the site and bond occupancy status for a site percolation
    process given a list of occupied sites.

    Parameters
    ----------
    ij : array_like
        An N x 2 array of [site_A, site_B] connections.  If two connected
        sites are both occupied they are part of the same cluster, as it
        the bond connecting them.

    occupied_sites : boolean
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
    occupied_bonds = np.all(occupied_sites[ij], axis=1)
    adj_mat = sprs.csr_matrix((occupied_bonds, (ij[:, 0], ij[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    clusters[~occupied_sites] = -1
    s_labels = spst.rankdata(clusters + 1, method="dense") - 1
    if np.any(~occupied_sites):
        s_labels -= 1
    b_labels = np.amin(s_labels[ij], axis=1)
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


def bond_percolation(ij, occupied_bonds):
    r"""
    Calculates the site and bond occupancy status for a bond percolation
    process given a list of occupied bonds.

    Parameters
    ----------
    ij : array_like
        An N x 2 array of [site_A, site_B] connections.  A site is
        considered occupied if any of it's connecting bonds are occupied.

    occupied_bonds: boolean
        A list indicating whether a bond is occupied or not

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to.  A value of -1 indicates uninvaded.

    Notes
    -----
    The ``connected_components`` function of scipy.sparse.csgraph will give ALL
    sites a cluster number whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    from collections import namedtuple
    Np = np.amax(ij) + 1
    # Find occupied sites based on occupied bonds
    # (the following 2 lines are not needed but worth keeping for future ref)
    # occupied_sites = np.zeros([Np, ], dtype=bool)
    # np.add.at(occupied_sites, ij[occupied_bonds].flatten(), True)
    adj_mat = sprs.csr_matrix((occupied_bonds, (ij[:, 0], ij[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    # Clusters of size 1 only occur if all a site's bonds are uninvaded
    valid_clusters = np.bincount(clusters) > 1
    mapping = -np.ones(shape=(clusters.max()+1, ), dtype=int)
    mapping[valid_clusters] = np.arange(0, valid_clusters.sum())
    s_labels = mapping[clusters]
    # Bond inherit the cluster number of its connected sites
    b_labels = np.amin(s_labels[ij], axis=1)
    # Set bond cluster to -1 if not actually occupied
    b_labels[~occupied_bonds] = -1
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


def find_clusters(network, mask=[], t_labels=False):
    r"""
    Identify connected clusters of pores in the network.  This method can
    also return a list of throat cluster numbers, which correspond to the
    cluster numbers of the pores to which the throat is connected.  Either
    site and bond percolation can be considered, see description of input
    arguments for details.

    Parameters
    ----------
    network : OpenPNM Network Object
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

    Examples
    --------
    >>> import openpnm as op
    >>> import numpy as np
    >>> pn = op.network.Cubic(shape=[25, 25, 1])
    >>> pn['pore.seed'] = np.random.rand(pn.Np)
    >>> pn['throat.seed'] = np.random.rand(pn.Nt)


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
    p_clusters = (clusters + 1)*(network.tomask(pores=Ps).astype(int)) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']][:, 0]
    # Label non-invaded throats with -1
    t_clusters[~tmask] = -1

    return (p_clusters, t_clusters)


def find_path(network, pore_pairs, weights=None):
    r"""
    Find the shortest path between pairs of pores.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object on which the search should be performed

    pore_pairs : array_like
        An N x 2 array containing N pairs of pores for which the shortest
        path is sought.

    weights : array_like, optional
        An Nt-long list of throat weights for the search.  Typically this
        would be the throat lengths, but could also be used to represent
        the phase configuration.  If no weights are given then the
        standard topological connections of the Network are used.

    Returns
    -------
    A dictionary containing both the pores and throats that define the
    shortest path connecting each pair of input pores.

    Notes
    -----
    The shortest path is found using Dijkstra's algorithm included in the
    scipy.sparse.csgraph module

    TODO: The returned throat path contains the correct values, but not
    necessarily in the true order

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[3, 3, 3])
    >>> a = op.topotools.find_path(network=pn, pore_pairs=[[0, 4], [0, 10]])
    >>> a['pores']
    [array([0, 1, 4]), array([ 0,  1, 10])]
    >>> a['throats']
    [array([ 0, 19]), array([ 0, 37])]
    """
    Ps = np.array(pore_pairs, ndmin=2)
    if weights is None:
        weights = np.ones_like(network.Ts)
    graph = network.create_adjacency_matrix(weights=weights, fmt='csr',
                                            drop_zeros=False)
    paths = csgraph.dijkstra(csgraph=graph, indices=Ps[:, 0],
                             return_predecessors=True)[1]
    pores = []
    throats = []
    for row in range(0, np.shape(Ps)[0]):
        j = Ps[row][1]
        ans = []
        while paths[row][j] > -9999:
            ans.append(j)
            j = paths[row][j]
        ans.append(Ps[row][0])
        ans.reverse()
        pores.append(np.array(ans, dtype=int))
        Ts = network.find_neighbor_throats(pores=ans, mode='xnor')
        throats.append(np.array(Ts, dtype=int))
    pdict = PrintableDict
    dict_ = pdict(**{'pores': pores, 'throats': throats})
    return dict_
