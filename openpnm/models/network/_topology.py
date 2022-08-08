r"""
Topology
--------
Pore-scale models related to topology of the network.

"""
from numpy.linalg import norm
from scipy.sparse import csgraph
import numpy as np

__all__ = [  # Keep this alphabetical for easier inspection of what's imported
    'coordination_number',
    'distance_to_furthest_neighbor',
    'distance_to_nearest_neighbor',
    'distance_to_nearest_pore',
    'pore_to_pore_distance',
    'reduce_coordination',
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


def reduce_coordination(network, z):
    r"""
    Deletes throats on network to match specified average coordination number

    Parameters
    ----------
    target : Network
        The network whose throats are to be trimmed
    z : scalar
        The desired average coordination number.  It is not possible to specify
        the distribution of the coordination, only the mean value.

    Returns
    -------
    trim : ndarray
        A boolean array with ``True`` values indicating which pores to trim
        (using ``op.topotools.trim``) to obtain the desired average
        coordination number.

    Notes
    -----
    This method first finds the minimum spanning tree of the network using
    random weights on each throat, then assures that these throats are *not*
    deleted, in order to maintain network connectivity.  The list of throats
    to trim is generated randomly from the throats *not* on the spanning tree.

    """
    # Find minimum spanning tree using random weights
    am = network.create_adjacency_matrix(weights=np.random.rand(network.Nt),
                                         triu=False)
    mst = csgraph.minimum_spanning_tree(am, overwrite=True)
    mst = mst.tocoo()

    # Label throats on spanning tree to avoid deleting them
    Ts = network.find_connecting_throat(mst.row, mst.col)
    Ts = np.hstack(Ts)
    network['throat.mst'] = False
    network['throat.mst'][Ts] = True

    # Trim throats not on the spanning tree to acheive desired coordination
    Ts = np.random.permutation(network.throats('mst', mode='nor'))
    del network['throat.mst']
    Ts = Ts[:int(network.Nt - network.Np*(z/2))]
    Ts = network.to_mask(throats=Ts)
    return Ts
