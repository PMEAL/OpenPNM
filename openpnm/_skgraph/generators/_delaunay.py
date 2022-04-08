import numpy as np
import scipy.spatial as sptl
from openpnm._skgraph import settings


def delaunay(points, shape=[1, 1, 1]):
    r"""
    Generate a network based on Delaunay triangulation of random points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate
    shape : array_like
        Indicates the size and shape of the domain

    Returns
    -------
    network : dict
        A dictionary containing 'vert.coords' and 'edge.conns'
    tri : Delaunay tessellation object
        The Delaunay tessellation object produced by ``scipy.spatial.Delaunay``
    """
    from openpnm.topotools import tri_to_am
    from openpnm.topotools.generators import tools

    node_prefix = settings.node_prefix
    edge_prefix = settings.edge_prefix

    points = tools.parse_points(points=points, shape=shape)
    mask = ~np.all(points == 0, axis=0)
    tri = sptl.Delaunay(points=points[:, mask])
    coo = tri_to_am(tri)
    d = {}
    d[node_prefix+'.coords'] = points
    d[edge_prefix+'.conns'] = np.vstack((coo.row, coo.col)).T
    return d, tri


if __name__ == "__main__":
    settings.node_prefix = 'node'
    settings.edge_prefix = 'edge'
    # Make a 2D network based on number of points
    dn, tri = delaunay(points=50, shape=[1, 1, 0])
    print(dn.keys())
    print(dn['node.coords'].shape)
    print(dn['edge.conns'].shape)
    # Make a 2D network based on number of points
    dn, tri = delaunay(points=50, shape=[1, 0, 1])
    print(dn.keys())
    print(dn['node.coords'].shape)
    print(dn['edge.conns'].shape)
    # Make a 3D network based on number of points
    dn, tri = delaunay(points=50, shape=[1, 1, 1])
    print(dn.keys())
    print(dn['node.coords'].shape)
    print(dn['edge.conns'].shape)
    # Make a 3D cylinder
    dn, tri = delaunay(points=50, shape=[1, 1])
    print(dn.keys())
    print(dn['node.coords'].shape)
    print(dn['edge.conns'].shape)
    # Make a 2D circle
    dn, tri = delaunay(points=50, shape=[1, 0])
    print(dn.keys())
    print(dn['node.coords'].shape)
    print(dn['edge.conns'].shape)
    # Make a 3D sphere
    dn, tri = delaunay(points=50, shape=[1])
    print(dn.keys())
    print(dn['node.coords'].shape)
    print(dn['edge.conns'].shape)
