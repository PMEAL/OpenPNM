"""
===============================================================================
Delaunay: Generate random networks based on the Delaunay Tessellation
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial as sptl
from openpnm.network import GenericNetwork
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Delaunay(GenericNetwork):
    r"""
    Generates a random network by performing at Delaunay Triangulation or
    tessellation on a given set of input points.

    Parameters
    ----------
    num_points : scalar
        The number of points to place in the domain, which will become the
        pore centers after the tessellation is performed.  This value is
        ignored if ``points`` are given.

    points : array_like
        An array of coordinates indicating the [x, y, z] locations of each
        point to use in the tessellation.

    Notes
    -----
    This class always performs the tessellation on the set of points, but
    trims any points that lie outside the given domain ``shape``.

    Examples
    --------
    >>> import openpnm as op
    >>> gn = op.network.Delaunay(num_points=50, shape=[1, 1, 0])


    """

    def __init__(self, num_points=None, shape=None, points=None, **kwargs):
        # Deal with input arguments
        if points is None:
            if num_points is None:
                raise Exception('Must specify either "points" or "num_points"')
            points = topotools.generate_base_points(num_points=num_points,
                                                    domain_size=shape)

        # Deal with points that are only 2D...they break Delaunay
        if points.shape[1] == 3 and len(sp.unique(points[:, 2])) == 1:
            points = points[:, :2]

        # Perform Delaunay tessellation on points
        tri = sptl.Delaunay(points)

        # Create an empty list-of-list matrix
        lil = sprs.lil_matrix((tri.npoints, tri.npoints))
        # Scan through Delaunay triangulation to retrieve pairs
        indices, indptr = tri.vertex_neighbor_vertices
        for k in range(tri.npoints):
            lil.rows[k] = indptr[indices[k]:indices[k+1]]
            lil.data[k] = sp.ones_like(lil.rows[k])

        # Convert to coo format
        coo = lil.tocoo()
        # Convert to normal conns
        conns = sp.vstack((coo.row, coo.col)).T
        # Convert to upper triangular
        conns = sp.sort(conns, axis=1)
        self['throat.conns'] = conns

        if points.shape[1] == 2:  # Make points 3D if necessary
            points = sp.vstack((points.T, sp.zeros((points.shape[0], )))).T
        self['pore.coords'] = points

        # coords and conns are needed to initialize GenericNetwork
        super().__init__(**kwargs)

        # Attach delaunay object in case someone needs it?
        self._tri = tri

        # Determine which throats constitute the Gabriel graph, for other uses
        # Find centroid of each pair of nodes
        c = points[conns]
        m = (c[:, 0, :] + c[:, 1, :])/2
        # Find radius of circle connecting each pair of nodes
        r = sp.sqrt(sp.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
        # Use KD-Tree to find distance to nearest neighbors
        tree = sptl.cKDTree(points)
        n = tree.query(x=m, k=1)[0]
        # Identify throats whose centroid is not near an unconnected node
        g = sp.around(n, decimals=5) == sp.around(r, decimals=5)
        # Label gabriel throats as such, for use elsewhere if needed
        self['throat.gabriel'] = g
