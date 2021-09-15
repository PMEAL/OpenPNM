import scipy.spatial as sptl
import scipy.sparse as sprs
from openpnm.network.generators import cubic
from openpnm.topotools import tri_to_am
import numpy as np
from numba import njit


@njit
def len_lil(lil):
    indptr = [len(i) for i in lil]
    return indptr


def bcc(shape, spacing, mode='kdtree'):
    shape = np.array(shape)
    spacing = np.array(spacing)
    net1 = cubic(shape=shape, spacing=1)
    net2 = cubic(shape=shape-1, spacing=1)
    net2['pore.coords'] += 0.5
    crds = np.concatenate((net1['pore.coords'], net2['pore.coords']))
    if mode.startswith('tri'):
        tri = sptl.Delaunay(points=crds)
        am = tri_to_am(tri)
        conns = np.vstack((am.row, am.col)).T
        # Trim diagonal connections between cubic pores
        L = np.sqrt(np.sum(np.diff(crds[conns], axis=1)**2, axis=2)).flatten()
        conns = conns[L <= 1]
    elif mode.startswith('kd'):
        tree1 = sptl.KDTree(crds)
        # Method 1
        hits = tree1.query_ball_point(crds, r=1)
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

    d = {}
    d['pore.coords'] = crds*spacing
    d['throat.conns'] = conns
    return d


if __name__ == '__main__':
    import openpnm as op
    pn = op.network.GenericNetwork()
    pn.update(bcc([3, 3, 3], 1, mode='tri'))
    pn['pore.all'] = np.ones((np.shape(pn.coords)[0]), dtype=bool)
    pn['throat.all'] = np.ones((np.shape(pn.conns)[0]), dtype=bool)
    fig = op.topotools.plot_connections(pn)
    fig = op.topotools.plot_coordinates(pn, ax=fig)
