from openpnm.topotools.generators import delaunay as _delaunay, tools
import scipy.spatial as sptl
import numpy as np


def gabriel(points=None, delaunay=None, shape=None):
    r"""
    Generate a network based on a Gabriel tessellation, which is a subset of
    the Delaunay triangulation

    Parameters
    ----------
    points : array_like or int, optional
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate.
        This can be omitted if ``delaunay`` is provided.
    delaunay : network dictionary, optional
        A dictionary containing 'pore.coords' and 'throat.conns' as produced
        by the ``delaunay`` function.  If ``points`` are provided this is
        ignored.
    shape : array_like
        Indicates the size and shape of the domain

    """
    if points is not None:
        delaunay, tri = _delaunay(points=points, shape=shape)
    # Find centroid or midpoint of each edge in conns
    c = delaunay['vert.coords'][delaunay['edge.conns']]
    m = (c[:, 0, :] + c[:, 1, :])/2
    # Find the radius sphere between each pair of nodes
    r = np.sqrt(np.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
    # Use the kd-tree function in Scipy's spatial module
    tree = sptl.cKDTree(delaunay['vert.coords'])
    # Find the nearest point for each midpoint
    n = tree.query(x=m, k=1)[0]
    # If nearest point to m is at distance r, then the edge is a Gabriel edge
    g = n >= r*(0.999)  # This factor avoids precision errors in the distances
    d = {}
    d.update(delaunay)
    # Reduce the connectivity to all True values found in g
    d['edge.conns'] = delaunay['edge.conns'][g]
    d['vert.coords'] = delaunay['vert.coords']
    return d


if __name__ == '__main__':
    from openpnm.topotools.generators import delaunay
    # Make a 2D network using points
    gb = gabriel(points=50, shape=[1, 1, 0])
    print(gb.keys())
    print(gb['vert.coords'].shape)
    print(gb['edge.conns'].shape)
    # Make a 2D network based on a pre-existing delaunay network
    dn, tri = delaunay(points=50, shape=[1, 1, 0])
    gb = gabriel(delaunay=dn)
    print(gb.keys())
    print(gb['vert.coords'].shape)
    print(gb['edge.conns'].shape)
    # Make a 3D network using points
    gb = gabriel(points=50, shape=[1, 1, 1])
    print(gb.keys())
    print(gb['vert.coords'].shape)
    print(gb['edge.conns'].shape)
    # Make a 3D network based on a pre-existing delaunay network
    dn, tri = delaunay(points=50, shape=[1, 1, 1])
    gb = gabriel(delaunay=dn)
    print(gb.keys())
    print(gb['vert.coords'].shape)
    print(gb['edge.conns'].shape)
