import numpy as np
import scipy.spatial as sptl
from openpnm._skgraph import settings
from openpnm._skgraph.tools import vor_to_am, isoutside
from openpnm._skgraph.generators import tools
from openpnm._skgraph.operations import trim_nodes


def voronoi(points, shape=[1, 1, 1], trim=True, tolerance=0.0,
            node_prefix='node', edge_prefix='edge'):
    r"""
    Generate a network based on a Voronoi tessellation of base points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used
        directly, or a scalar indicating the number of points to generate.
    shape : array_like
        The size and shape of the domain
    trim : boolean
        If ``True`` (default) then any vertices laying outside the domain
        given by ``shape`` are removed (as are the edges connected to them).
    tolerance : float
        The tolerance to use when deciding if a point is within the domain or
        not when ``trim=True``. The default is 0.  Increasing this is usually
        only needed when dealing with curved domain faces (spheres and
        cylinders).

    Returns
    -------
    network : dict
        A dictionary containing node coordinates and edge connections
    vor : Voronoi tessellation object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``

    """
    points = tools.parse_points(points=points, shape=shape)
    mask = ~np.all(points == 0, axis=0)
    # Perform tessellation
    vor = sptl.Voronoi(points=points[:, mask])
    # Convert to adjecency matrix
    coo = vor_to_am(vor)
    # Write values to dictionary
    d = {}
    conns = np.vstack((coo.row, coo.col)).T
    d[edge_prefix+'.conns'] = conns

    # Convert coords to 3D if necessary
    # Rounding is crucial since some voronoi verts endup outside domain
    pts = np.around(vor.vertices, decimals=10)
    if mask.sum() < 3:
        coords = np.zeros([pts.shape[0], 3], dtype=float)
        coords[:, mask] = pts
    else:
        coords = pts

    d[node_prefix+'.coords'] = coords

    if trim:
        hits = isoutside(d, shape=shape, rtol=tolerance)
        d = trim_nodes(d, hits)

    return d, vor


if __name__ == "__main__":
    from openpnm._skgraph.visualization import plot_edges
    vn, vor = voronoi(points=50, shape=[1, 0, 1])
    print(vn.keys())
    print(vn['node.coords'].shape)
    print(vn['edge.conns'].shape)

    shape = [1, 0]
    pts = tools.parse_points(points=1000, shape=shape, reflect=True)
    vn, vor = voronoi(points=pts, shape=shape, trim=True, tolerance=[0.2, 0])
    plot_edges(vn)
