from numpy.linalg import norm as _norm
import numpy as _np


r"""
Pore-scale models related to topology of the network.
"""


def coordination_number(target):
    r"""
    Find the number of neighbors for each pore
    """
    network = target.network
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    vals = N[network.pores(target.name)]
    return vals


def pore_to_pore_distance(target):
    r"""
    Find the center to center distance between each pair of pores
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = _norm(C1 - C2, axis=1)
    return values


def distance_to_nearest_neighbor(target):
    r"""
    Find the distance between each pore and its closest topological neighbor
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = _norm(C1 - C2, axis=1)

    data = values
    im = network.create_incidence_matrix()
    values = _np.ones((network.Np, ))*_np.inf
    _np.minimum.at(values, im.row, data[im.col])
    return _np.array(values)


def distance_to_furthest_neighbor(target):
    r"""
    Find the distance between each pore and its furthest topological neighbor
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = _norm(C1 - C2, axis=1)

    data = values
    im = network.create_incidence_matrix()
    values = _np.zeros((network.Np, ))
    _np.maximum.at(values, im.row, data[im.col])
    return _np.array(values)


def cluster_number(target):
    r"""
    Assign a cluster number to each pore
    """
    net = target.network
    from scipy.sparse import csgraph as csg
    am = net.create_adjacency_matrix(fmt='coo', triu=True)
    N, Cs = csg.connected_components(am, directed=False)
    return Cs


def cluster_size(target, cluster=None):
    r"""
    Find the size of the cluster to which each pore belongs

    Parameters
    ----------
    network : dict
        The OpenPNM network object
    cluster : str, optional
        Dict key pointing to the array containing the cluster number of each
        pore.  If not provided then it will be calculated.

    Returns
    -------
    cluster_size : ndarray
        An Np-long array containing the size of the cluster to which each pore
        belongs

    """
    net = target.network
    if cluster is None:
        from scipy.sparse import csgraph as csg
        am = net.create_adjacency_matrix(fmt='coo', triu=True)
        N, cluster_num = csg.connected_components(am, directed=False)
    else:
        cluster_num = net[cluster]
    Cs, ind, N = _np.unique(cluster_num, return_inverse=True, return_counts=True)
    values = N[ind]
    return values


def isolated_pores(target):
    r"""
    find which pores, if any, are not connected to a throat
    """
    net = target.network
    values = _np.ones(net.Np, dtype=bool)
    hits = _np.unique(net.conns)
    values[hits] = False
    return values


def reversed_throats(target):
    r"""
    Find any throat connections that are pointing from j -> i where j > i
    """
    net = target.network
    hits = net.conns[:, 0] > net.conns[:, 1]
    return hits


def looped_throats(target):
    r"""
    Find any throats that are connected to the same pore on both ends
    """
    net = target.network
    hits = net.conns[:, 0] == net.conns[:, 1]
    return hits


def headless_throats(target):
    r"""
    Find any throats that point to a non-existent pore
    """
    net = target.network
    hits = _np.any(net.conns > (net.Np -1), axis=1)
    return hits


def duplicate_throats(target):
    r"""
    Find repeat occurrences of throat connections
    """
    net = target.network
    conns = net.conns
    iconns = conns[:, 0] + 1j*conns[:, 1]
    hits, inds = _np.unique(iconns, return_inverse=True)
    values = _np.ones(net.Nt, dtype=bool)
    values[inds] = False
    return values


def distance_to_nearest_pore(target):
    r"""
    Find distance to and index of nearest pore even if not topologically
    connected
    """
    import scipy.spatial as sptl
    net = target.network
    coords = net.coords
    tree = sptl.KDTree(coords)
    ds, ids = tree.query(coords, k=2)
    values = ds[:, 1]
    return values


def count_coincident_pores(target, thresh=1e-6):
    r"""
    Count number of pores that are spatially coincident with other pores

    Parameters
    ----------
    network : dict
        The OpenPNM network object
    thresh : float
        The distance below which two pores are considered spatially coincident

    Returns
    -------
    count : ndarray
        A numpy array of Np length containing the number coincident pores
    """
    # This needs to be a bit complicated because it cannot be assumed
    # the coincident pores are topologically connected
    import scipy.spatial as sptl
    net = target.network
    coords = net.coords
    tree = sptl.KDTree(coords)
    hits = tree.query_pairs(r=thresh)
    arr = _np.array(list(hits)).flatten()
    v, n = _np.unique(arr, return_counts=True)
    values = _np.zeros(net.Np, dtype=int)
    values[v.astype(int)] = n
    return values


def find_coincident_pores(target, thresh=1e-6):
    r"""
    Find the indices of coincident pores

    Parameters
    ----------
    network : dict
        The OpenPNM network object
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
    network = target.network
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
