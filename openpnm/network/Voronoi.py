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
    By definition the given points will each lie in the center of a Voronoi
    cell, so they will not be the pore centers. The pores in the returned
    network are actually the vertices of the Voronoi tessellation, thus will
    differ from the number of points supplied.

    """

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        # Clean-up input points
        super().__init__(**kwargs)
        points = self._parse_points(shape=shape, points=points)
        network, vor = voronoi(points=points, shape=shape)
        # Initialize network object
        self.update(network)
        topotools.trim(network=self, pores=self.pores('crop'))
        del self['pore.crop']
