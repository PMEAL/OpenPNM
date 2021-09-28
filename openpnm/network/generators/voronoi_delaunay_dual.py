import scipy.spatial as sptl
import scipy.sparse as sprs
from openpnm import topotools
from openpnm.network.generators import tools
import numpy as np


def voronoi_delaunay_dual(points, shape):
    r"""
    Generate a dual Voronoi-Delaunay network from given base points

    Parameters
    ----------
    points : array_like or scalar
        The points to be tessellated.  If a scalar is given a set of points
        of that size is generated.
    shape : array_like
        The size of the domain in which the points lie.

    Returns
    -------
    network : dict
        A dictionary containing ``'pore.coords'`` and ``'throat.conns'``.
    vor : Voronoi object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``
    tri : Delaunay object
        The Delaunay triangulation object produced ``scipy.spatial.Delaunay``

    """
    points = tools.parse_points(points=points, shape=shape)

    # Perform tessellation
    vor = sptl.Voronoi(points=points)
    tri = sptl.Delaunay(points=points)

    # Combine points
    pts_all = np.vstack((vor.points, vor.vertices))
    Nall = np.shape(pts_all)[0]

    # Create adjacency matrix in lil format for quick construction
    am = sprs.lil_matrix((Nall, Nall))
    for ridge in vor.ridge_dict.keys():
        # Make Delaunay-to-Delauny connections
        for i in ridge:
            am.rows[i].extend([ridge[0], ridge[1]])
        # Get voronoi vertices for current ridge
        row = vor.ridge_dict[ridge].copy()
        # Index Voronoi vertex numbers by number of delaunay points
        row = [i + vor.npoints for i in row if i > -1]
        # Make Voronoi-to-Delaunay connections
        for i in ridge:
            am.rows[i].extend(row)
        # Make Voronoi-to-Voronoi connections
        row.append(row[0])
        for i in range(len(row)-1):
            am.rows[row[i]].append(row[i+1])

    # Finalize adjacency matrix by assigning data values
    am.data = am.rows  # Values don't matter, only shape, so use 'rows'
    # Convert to COO format for direct acces to row and col
    am = am.tocoo()
    # Extract rows and cols
    conns = np.vstack((am.row, am.col)).T

    # Convert to sanitized adjacency matrix
    am = topotools.conns_to_am(conns)
    # Finally, retrieve conns back from am
    conns = np.vstack((am.row, am.col)).T

    # Translate adjacency matrix and points to OpenPNM format
    coords = np.around(pts_all, decimals=10)
    if coords.shape[1] == 2:  # Make points back into 3D if necessary
        coords = np.vstack((coords.T, np.zeros((coords.shape[0], )))).T

    network = {}
    network['pore.coords'] = coords
    network['throat.conns'] = conns

    return network, vor, tri