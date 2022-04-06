import scipy.spatial as sptl
import scipy.sparse as sprs
from openpnm._skgraph.generators import cubic
import numpy as np
from numba import njit


@njit
def len_lil(lil):
    indptr = [len(i) for i in lil]
    return indptr


def fcc(shape, spacing=1, mode='kdtree', node_prefix='node', edge_prefix='edge'):
    r"""
    Generate a face-centered cubic lattice

    Parameters
    ----------
    shape : array_like
        The number of unit cells in each direction.  A unit cell has vertices
        on all 8 corners, 6 faces and a single vertex at its center.
    spacing : array_like or float
        The size of a unit cell in each direction. If an scalar is given it is
        applied in all 3 directions.
    mode : str
        Dictate how neighbors are found.  Options are:

        ===============  =====================================================
        mode             meaning
        ===============  =====================================================
        'kdtree'         Uses ``scipy.spatial.KDTree`` to find all neighbors
                         within the unit cell.
        'triangulation'  Uses ``scipy.spatial.Delaunay`` to find all neighbors.
        ===============  =====================================================

    node_prefix : str, optional
        If a custom prefix is used to indicate node arrays, such as ('site', or
        'vertex') it can be specified here.  The defaul it 'node'.
    edge_prefix : str, optional
        If a custom prefix is used to indicate site arrays, such as ('bond', or
        'link') it can be specified here.  The defaul it 'edge'.


    Returns
    -------
    network : dict
        A dictionary containing 'coords' and 'conns'

    Notes
    -----
    It is not clear whether KDTree of Delaunay are faster. In fact it is
    surely possible to find the neighbors formulaically but this is not
    implemented yet.

    """
    from openpnm.topotools import tri_to_am
    shape = np.array(shape)
    # Create base cubic network of corner sites
    net1 = cubic(shape=shape)
    # Create 3 networks to become face sites
    net2 = cubic(shape=shape - [1, 1, 0])
    net3 = cubic(shape=shape - [1, 0, 1])
    net4 = cubic(shape=shape - [0, 1, 1])
    # Offset pore coords by 1/2 a unit cell
    net2['node.coords'] += np.array([0.5, 0.5, 0])
    net3['node.coords'] += np.array([0.5, 0, 0.5])
    net4['node.coords'] += np.array([0, 0.5, 0.5])
    crds = np.concatenate((net1['node.coords'],
                           net2['node.coords'],
                           net3['node.coords'],
                           net4['node.coords']))
    corner_labels = np.concatenate(
        (np.ones(net1['node.coords'].shape[0], dtype=bool),
         np.zeros(net2['node.coords'].shape[0], dtype=bool),
         np.zeros(net3['node.coords'].shape[0], dtype=bool),
         np.zeros(net4['node.coords'].shape[0], dtype=bool)))
    if mode.startswith('tri'):
        tri = sptl.Delaunay(points=crds)
        am = tri_to_am(tri)
        conns = np.vstack((am.row, am.col)).T
        # Trim diagonal connections between cubic pores
        L = np.sqrt(np.sum(np.diff(crds[conns], axis=1)**2, axis=2)).flatten()
        conns = conns[L <= 0.75]
    elif mode.startswith('kd'):
        tree1 = sptl.KDTree(crds)
        # Method 1
        hits = tree1.query_ball_point(crds, r=.75)
        # Method 2: Not sure which is faster
        # tree2 = sptl.KDTree(crds)
        # hits = tree1.query_ball_tree(tree1, r=1)
        indices = np.hstack(hits)
        # Convert to CSR matrix
        indptr = [len(i) for i in hits]
        indptr.insert(0, 0)
        indptr = np.cumsum(indptr)
        am = sprs.csr_matrix((np.ones_like(indices), indices, indptr))
        am = sprs.triu(am, k=1)
        am = am.tocoo()
        conns = np.vstack((am.row, am.col)).T
    conns = np.vstack((net1['edge.conns'], conns))

    d = {}
    d[node_prefix+'.corner'] = corner_labels
    d[node_prefix+'.face'] = ~corner_labels
    d[node_prefix+'.coords'] = crds*spacing
    d[edge_prefix+'.conns'] = conns
    return d


if __name__ == '__main__':
    import openpnm as op
    import matplotlib.pyplot as plt
    net = fcc([3, 3, 3], 1, mode='tri')
    net['pore.coords'] = net.pop('node.coords')
    net['throat.conns'] = net.pop('edge.conns')
    net['pore.corner'] = net.pop('node.corner')
    net['pore.face'] = net.pop('node.face')
    pn = op.network.GenericNetwork()
    pn.update(net)
    pn['pore.all'] = np.ones((np.shape(pn.coords)[0]), dtype=bool)
    pn['throat.all'] = np.ones((np.shape(pn.conns)[0]), dtype=bool)
    fig, ax = plt.subplots()
    op.topotools.plot_connections(pn, ax=ax)
    op.topotools.plot_coordinates(pn, ax=ax)
