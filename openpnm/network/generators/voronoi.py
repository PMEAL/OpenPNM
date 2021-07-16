import scipy.spatial as sptl
import numpy as np
from openpnm.topotools import vor_to_am, isoutside


def voronoi(points, shape=[1, 1, 1], crop=True):
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
    return d
