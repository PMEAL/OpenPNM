r"""
Topology
--------
Pore-scale models related to topology of the network.

"""
from numpy.linalg import norm
import numpy as np

__all__ = [  # Keep this alphabetical for easier inspection of what's imported
    'coordination_number',
    'distance_to_furthest_neighbor',
    'distance_to_nearest_neighbor',
    'distance_to_nearest_pore',
    'pore_to_pore_distance',
]


def coordination_number(network):
    r"""
    Find the number of neighbors for each pore
    """
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    return N


def pore_to_pore_distance(network):
    r"""
    Find the center to center distance between each pair of pores
    """
    cn = network['throat.conns']
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = norm(C1 - C2, axis=1)
    return values


def distance_to_nearest_neighbor(network):
    r"""
    Find the distance between each pore and its closest topological neighbor
    """
    cn = network['throat.conns']
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    D = norm(C1 - C2, axis=1)
    im = network.create_incidence_matrix()
    values = np.ones((network.Np, ))*np.inf
    np.minimum.at(values, im.row, D[im.col])
    return np.array(values)


def distance_to_furthest_neighbor(network):
    r"""
    Find the distance between each pore and its furthest topological neighbor
    """
    throats = network.throats(network.name)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    D = norm(C1 - C2, axis=1)
    im = network.create_incidence_matrix()
    values = np.zeros((network.Np, ))
    np.maximum.at(values, im.row, D[im.col])
    return np.array(values)


def distance_to_nearest_pore(network):
    r"""
    Find distance to and index of nearest pore even if not topologically
    connected
    """
    import scipy.spatial as sptl
    coords = network.coords
    tree = sptl.KDTree(coords)
    ds, ids = tree.query(coords, k=2)
    values = ds[:, 1]
    return values
