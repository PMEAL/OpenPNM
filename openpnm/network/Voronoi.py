"""
===============================================================================
Voronoi: Generate random networks based on the Voronoi Tessellation
===============================================================================

"""

from openpnm.network import DelaunayVoronoiDual
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Voronoi(DelaunayVoronoiDual):
    r"""

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
