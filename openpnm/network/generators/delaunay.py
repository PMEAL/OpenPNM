import numpy as np
import scipy.spatial as sptl
import scipy.sparse as sprs


def delaunay(points, shape=[1, 1, 1]):
    shape = np.array(shape)
    if isinstance(points, int):
        points = np.random.rand(points, len(shape))*shape
    tri = sptl.Delaunay(points=points)
    indices, indptr = tri.vertex_neighbor_vertices
    coo = [[],[]]
    for k in range(tri.npoints):
        col = indptr[indices[k]:indices[k+1]]
        coo[0].extend(np.ones_like(col)*k)
        coo[1].extend(col)
    coo = sprs.coo_matrix((np.ones_like(coo[0]), (coo[0], coo[1])))
    coo = sprs.triu(coo, k=1)
    d = {}
    d['pore.coords'] = points
    d['throat.conns'] = np.vstack((coo.row, coo.col)).T
    return d
