"""
===============================================================================
Delaunay: Generate random networks based on the Delaunay Tessellation
===============================================================================

"""
from openpnm.network import DelaunayVoronoiDual
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Delaunay(DelaunayVoronoiDual):
    r"""
    Generates a random network by performing at Delaunay triangulation (aka
    tessellation) on a set of base points.

    Parameters
    ----------
    num_points : scalar
        The number of points to place in the domain, which will become the
        pore centers after the tessellation is performed.  This value is
        ignored if ``points`` are given.

    points : array_like
        An array of coordinates indicating the [x, y, z] locations of each
        point to use in the tessellation.  Note that the points must be given
        in rectilinear coordinates regardless of which domain ``shape`` was
        specified.  To convert between coordinate systems see the
        ``convert_coords`` function in the ``openpnm.topotools`` module.

    shape : array_like
        The size of the domain.  It's possible to create cubic, cylindrical,
        or spherical domains, as well as 2D square and circular by changing
        the domain ``shape`` as follows:

        [x, y, z] - will produce a normal cubic domain of dimension x, and
        and z

        [x, y, 0] - will produce a 2D square domain of size x by y

        [r, z] - will produce a cylindrical domain with a radius of r and
        height of z

        [r, 0] - will produce a 2D circular domain with a radius of r

        [r] - will produce a spherical domain with a radius of r

    Notes
    -----
    This class always performs the tessellation on the full set of points, then
    trims any points that lie outside the given domain ``shape``.  This is
    important for cases where base points have been reflected about the domain
    edges since all reflected points are deleted to reveal the smoothly
    tessellated surface.

    Examples
    --------
    >>> import openpnm as op
    >>> gn = op.network.Delaunay(num_points=50, shape=[1, 1, 0])

    """

    def __init__(self, shape=None, num_points=None, points=None, **kwargs):
        # Clean-up input points
        points = self._parse_points(shape=shape,
                                    num_points=num_points,
                                    points=points)
        # Initialize network object
        super().__init__(shape=shape, points=points, **kwargs)
        topotools.trim(network=self, pores=self.pores(['voronoi']))
        pop = ['pore.voronoi', 'throat.voronoi', 'throat.interconnect',
               'pore.delaunay', 'throat.delaunay']
        for item in pop:
            del self[item]
