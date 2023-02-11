r"""
Topology
--------
Pore-scale models related to topology of the network.

"""
from numpy.linalg import norm
import numpy as np
import scipy.spatial as sptl


__all__ = [  # Keep this alphabetical for easier inspection of what's imported
    'coordination_number',
    'distance_to_furthest_neighbor',
    'distance_to_nearest_neighbor',
    'distance_to_nearest_pore',
    'gabriel_edges',
    'pore_to_pore_distance',
]


def gabriel_edges(network):
    r"""
    Find throats which make a Gabriel subgraph

    Returns
    -------
    throats : ndarray
        An ndarray of boolean values with ``True`` indicating that a throat
        satisfies the conditions of Gabriel graph, meaning that a circle (or
        sphere) can be drawn between its two connected pores that does not
        contain any other pores.

    Notes
    -----
    Technically this should only be used on a Delaunay network, but it will
    work on any graph. By deleting all throats that are *not* identified by
    this fuction one would obtain the Gabriel graph [1].

    References
    ----------
    [1] `Wikipedia <https://en.wikipedia.org/wiki/Gabriel_graph>`_

    """
    dn = network
    # Find centroid or midpoint of each edge in conns
    c = dn['pore.coords'][dn['throat.conns']]
    m = (c[:, 0, :] + c[:, 1, :])/2
    # Find the radius the sphere between each pair of nodes
    r = np.sqrt(np.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
    # Use the kd-tree function in Scipy's spatial module
    tree = sptl.cKDTree(dn['pore.coords'])
    # Find the nearest point for each midpoint
    n = tree.query(x=m, k=1)[0]
    # If nearest point to m is at distance r, then the edge is a Gabriel edge
    g = n >= r*(0.999)  # This factor avoids precision errors in the distances
    return g


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
