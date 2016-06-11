"""
===============================================================================
Voronoi: Generate a random network based on Voronoi tessellation of random
points
===============================================================================

"""

import scipy as sp
import scipy.spatial as sptl
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Voronoi(GenericNetwork):
    r"""

    """

    def __init__(self, num_cells, domain_size=[1, 1, 1], base_points=None,
                 flat_edges=True, **kwargs):
        super().__init__(**kwargs)
        domain_size = sp.array(domain_size, ndmin=1)
        if base_points is None:
            base_points = sp.rand(num_cells, 3)
            base_points = base_points*domain_size
        self.vor = sptl.Voronoi(points=base_points)
