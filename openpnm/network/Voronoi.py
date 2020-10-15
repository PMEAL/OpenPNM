from openpnm.network import DelaunayVoronoiDual
from openpnm import topotools
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Voronoi(DelaunayVoronoiDual):
    r"""
    Random network formed by Voronoi tessellation of arbitrary base points

    Parameters
    ----------
    points : array_like, optional
        The base points around which to generate the Voronoi tessellation.

    num_points : scalar, optional
        If ``points`` is not supplied, then this must be given.  A sent of
        randomly located points will be generated.

    shape : array_like
        The size of the domain.  It's possible to create cubic as well as 2D
        square domains by changing the ``shape`` as follows:

        [x, y, z] - will produce a normal cubic domain of dimension x, and
        and z

        [x, y, 0] - will produce a 2D square domain of size x by y

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.

    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a *Project*.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a *Project* use ``openpnm.Project()``.

    Notes
    -----
    By definition these points will each lie in the center of a Voronoi cell,
    so they will not be the pore centers.  The number of pores in the
    returned network thus will differ from the number of points supplied

    """

    def __init__(self, shape=[1, 1, 1], num_points=None, points=None, **kwargs):
        # Clean-up input points
        points = self._parse_points(shape=shape,
                                    num_points=num_points,
                                    points=points)
        super().__init__(shape=shape, points=points, **kwargs)
        # Initialize network object
        topotools.trim(network=self, pores=self.pores('delaunay'))
        pop = ['pore.delaunay', 'throat.delaunay', 'throat.interconnect']
        for item in pop:
            del self[item]
