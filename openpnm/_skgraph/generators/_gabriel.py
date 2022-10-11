import numpy as np
import scipy.spatial as sptl
from openpnm._skgraph.generators import delaunay as _delaunay


__all__ = [
    'gabriel',
]


def gabriel(points=None, delaunay=None, shape=None,
            node_prefix='node', edge_prefix='edge'):
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
        A dictionary containing coords and conns as produced by the
        ``delaunay`` function.  If ``points`` are provided this is
        ignored.
    shape : array_like
        Indicates the size and shape of the domain

    Returns
    -------
    network : dict
        A dictionary containing 'node.coords' and 'edge.conns'

    """
    if points is not None:
        dn, tri = _delaunay(points=points, shape=shape,
                            node_prefix=node_prefix, edge_prefix=edge_prefix)
    else:
        dn = delaunay
    # Find centroid or midpoint of each edge in conns
    c = dn[node_prefix+'.coords'][dn[edge_prefix+'.conns']]
    m = (c[:, 0, :] + c[:, 1, :])/2
    # Find the radius the sphere between each pair of nodes
    r = np.sqrt(np.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
    # Use the kd-tree function in Scipy's spatial module
    tree = sptl.cKDTree(dn[node_prefix+'.coords'])
    # Find the nearest point for each midpoint
    n = tree.query(x=m, k=1)[0]
    # If nearest point to m is at distance r, then the edge is a Gabriel edge
    g = n >= r*(0.999)  # This factor avoids precision errors in the distances
    d = {}
    d.update(dn)
    # Reduce the connectivity to all True values found in g
    d[edge_prefix+'.conns'] = dn[edge_prefix+'.conns'][g]
    d[node_prefix+'.coords'] = dn[node_prefix+'.coords']
    return d


if __name__ == '__main__':
    from openpnm._skgraph.generators import delaunay
    # Make a 2D network using points
    gb = gabriel(points=50, shape=[1, 1, 0])
    print(gb.keys())
    print(gb['node.coords'].shape)
    print(gb['edge.conns'].shape)
    # Make a 2D network based on a pre-existing delaunay network
    dn, tri = delaunay(points=50, shape=[1, 1, 0])
    gb = gabriel(delaunay=dn)
    print(gb.keys())
    print(gb['node.coords'].shape)
    print(gb['edge.conns'].shape)
    # Make a 3D network using points
    gb = gabriel(points=50, shape=[1, 1, 1])
    print(gb.keys())
    print(gb['node.coords'].shape)
    print(gb['edge.conns'].shape)
    # Make a 3D network based on a pre-existing delaunay network
    dn, tri = delaunay(points=50, shape=[1, 1, 1])
    gb = gabriel(delaunay=dn)
    print(gb.keys())
    print(gb['node.coords'].shape)
    print(gb['edge.conns'].shape)
