from openpnm.network import DelaunayVoronoiDual
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Voronoi(DelaunayVoronoiDual):
    r"""
    Creates a network based on a Voronoi tessellation of arbitrary base points.

    Parameters
    ----------
    points : array_like, optional
        The base points around which to generate the Voronoi tessellation.

    num_points : scalar, optional
        If ``points`` is not supplied, then this must be given.  A sent of
        randomly located points will be generated.

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
    By definition these points will each lie in the center of a Voronoi cell,
    so they will not be the pore centers.  The number of pores in the
    returned network thus will differ from the number of points supplied

    """
    def __init__(self, shape=None, num_points=None, points=None, **kwargs):
        # Clean-up input points
        points = self._parse_points(shape=shape,
                                    num_points=num_points,
                                    points=points)
        # Initialize network object
        super().__init__(shape=shape, points=points, **kwargs)
        topotools.trim(network=self, pores=self.pores('delaunay'))
        pop = ['pore.delaunay', 'throat.delaunay', 'throat.interconnect']
        for item in pop:
            del self[item]
