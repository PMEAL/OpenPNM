import numpy as np
import scipy.spatial as sptl
from openpnm._skgraph import settings


def voronoi(points, shape=[1, 1, 1]):
    r"""
    Generate a network based on a Voronoi tessellation of base points

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
        A dictionary containing 'node.coords' and 'edge.conns'
    vor : Voronoi tessellation object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``

    """
    from openpnm.topotools import vor_to_am, isoutside
    from openpnm.topotools.generators import tools

    node_prefix = settings.node_prefix
    edge_prefix = settings.edge_prefix

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
    if np.any(mask == False):
        verts = np.zeros([np.shape(vor.vertices)[0], 3])
        for i, col in enumerate(np.where(mask)[0]):
            verts[:, col] = vor.vertices[:, i]
    else:
        verts = np.copy(vor.vertices)
    d[node_prefix+'.coords'] = verts

    return d, vor


if __name__ == "__main__":
    settings.node_prefix = 'node'
    settings.edge_prefix = 'edge'

    vn, vor = voronoi(points=50, shape=[1, 0, 1])
    print(vn.keys())
    print(vn['node.coords'].shape)
    print(vn['edge.conns'].shape)
    vn, vor = voronoi(points=50, shape=[1, 0, 1])
    print(vn.keys())
    print(vn['node.coords'].shape)
    print(vn['edge.conns'].shape)
