import numpy as np
import scipy.spatial as sptl
from openpnm.topotools import tri_to_am


def delaunay(points, shape=[1, 1, 1], ):
    r"""
    Generate a network based on Delaunay triangulation of random points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate.
    shape : array_like
        Indicates the size and shape of the domain.

    Returns
    -------
    network : dict
        A dictionary containing 'pore.coords' and 'throat.conns'
    tri : Delaunay tessellation object
        The Delaunay tessellation object produced by ``scipy.spatial.Delaunay``
    """
    shape = np.array(shape)
    if isinstance(points, int):
        points = np.random.rand(points, len(shape))*shape
    tri = sptl.Delaunay(points=points)
    coo = tri_to_am(tri)
    d = {}
    d['pore.coords'] = points
    d['throat.conns'] = np.vstack((coo.row, coo.col)).T
    return d, tri
