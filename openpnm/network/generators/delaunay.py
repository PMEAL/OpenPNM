import numpy as np
import scipy.spatial as sptl
from openpnm.topotools import tri_to_am


def delaunay(points, shape=[1, 1, 1]):
    shape = np.array(shape)
    if isinstance(points, int):
        points = np.random.rand(points, len(shape))*shape
    tri = sptl.Delaunay(points=points)
    coo = tri_to_am(tri)
    d = {}
    d['pore.coords'] = points
    d['throat.conns'] = np.vstack((coo.row, coo.col)).T
    return d
