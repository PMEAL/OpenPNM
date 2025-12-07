import numpy as np
import scipy.spatial as sptl
from openpnm._skgraph.tools import vor_to_am, isoutside
from openpnm._skgraph.generators import tools
from openpnm._skgraph.operations import trim_nodes


def voronoi(
    points,
    shape=[1, 1, 1],
    trim=True,
    reflect=False,
    f=1,
    relaxation=0,
    node_prefix='node',
    edge_prefix='edge',
):
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
    reflect : boolean, optional (default = ``False``)
        If ``True`` then points are reflected across each face of the domain
        prior to performing the tessellation. Enabling this behavior creates
        flat faces on all sizes of the domain.
    f : float
        The fraction of points which should be reflected.  The default is 1 which
        reflects all the points in the domain, but this can lead to a lot of
        unnecessary points, so setting to 0.1 or 0.2 helps speed, but risks that
        the tessellation may not have smooth faces if not enough points are
        reflected.
    relaxation : int, optional (default = 0)
        The number of iterations to use for relaxing the base points. This is
        sometimes called `Lloyd's algorithm
        <https://en.wikipedia.org/wiki/Lloyd%27s_algorithm>`_. This function computes
        the new base points as the simple average of the Voronoi vertices instead
        of rigorously finding the center of mass, which is quite time consuming.
        To use the rigorous method, call the ``lloyd_relaxation`` function manually
        to obtain relaxed points, then pass the points directly to this funcion.
        The results are quite stable after only a few iterations.

    Returns
    -------
    network : dict
        A dictionary containing node coordinates and edge connections
    vor : Voronoi tessellation object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``

    """
    points = tools.parse_points(points=points, shape=shape, reflect=reflect, f=f)
    mask = ~np.all(points == 0, axis=0)
    # Perform tessellation
    vor = sptl.Voronoi(points=points[:, mask])
    for _ in range(relaxation):
        points = tools.lloyd_relaxation(vor, mode='rigorous')
        # Reparse points
        d = {}
        d[node_prefix+'.coords'] = points
        keep = ~isoutside(network=d, shape=shape)
        points = points[keep]
        points = tools.parse_points(points=points, shape=shape, reflect=reflect)
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
        hits = isoutside(d, shape=shape)
        d = trim_nodes(d, hits)
    return d, vor
