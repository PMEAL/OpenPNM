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
    def __init__(self, num_points=None, shape=None, points=None,
                 trim_domain=False, **kwargs):

        super().__init__(num_points=num_points, shape=shape, points=points,
                         trim_domain=trim_domain, **kwargs)
        topotools.trim(self, pores=self.pores('delaunay'))
        del self['throat.interconnect']
        del self['throat.delaunay']
        del self['pore.delaunay']
        del self['pore.voronoi']
        del self['throat.voronoi']
