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


def nearest_neighbor_distance(target):
    r"""
    Find the distance between each pore and its closest direct neighbor
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


def furthest_neighbor_distance(target):
    r"""
    Find the distance between each pore and its furthest direct neighbor
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
    from scipy.sparse import csgraph as csg
    net = target.network
    am = net.create_adjacency_matrix(fmt='coo', triu=True)
    N, Cs = csg.connected_components(am, directed=False)[1]
    return Cs


def cluster_size(target, cluster='pore.cluster_number'):
    r"""
    Find the size of the cluster to which each pore belongs
    """
    net = target.network
    Cs, ind, N = _np.unique(net[cluster], return_inverse=True, return_counts=True)
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
    Finds repeat occurrences of throat connections
    """
    net = target.network
    conns = net.conns
    iconns = conns[:, 0] + 1j*conns[:, 1]
    hits, inds = _np.unique(iconns, return_inverse=True)
    values = _np.ones(net.Nt, dtype=bool)
    values[inds] = False









