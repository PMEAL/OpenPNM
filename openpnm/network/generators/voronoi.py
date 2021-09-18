import scipy.spatial as sptl
import numpy as np
from openpnm.topotools import vor_to_am, isoutside


def voronoi(points, shape=[1, 1, 1], crop=True):
    r"""
    Generate a network based on a Voronoi tessellation of random points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate.
    shape : array_like
        Indicates the size and shape of the domain.
    crop : boolean
        If ``True`` (default) any points in the tessellation that lie outside
        the domain (defined by ``shape``) are removed.

    Returns
    -------
    network : dict
        A dictionary containing 'pore.coords' and 'throat.conns'
    vor : Voronoi tessellation object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``

    """
    shape = np.array(shape)
    if isinstance(points, int):
        points = np.random.rand(points, len(shape))*shape

    # Perform tessellation
    vor = sptl.Voronoi(points=points)
    # Convert to adjecency matrix
    coo = vor_to_am(vor)
    # Write values to dictionary
    d = {}
    conns = np.vstack((coo.row, coo.col)).T
    d['throat.conns'] = conns
    d['pore.coords'] = vor.vertices
    if crop:
        d['pore.crop'] = isoutside(vor.vertices, shape=shape)
    return d, vor
