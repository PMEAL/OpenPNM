"""
===============================================================================
Voronoi: Generate random networks based on the Voronoi Tessellation
===============================================================================

"""

import scipy as sp
import scipy.spatial as sptl
from openpnm.network import GenericTessellation
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Voronoi(GenericTessellation):
    r"""

    """
    def __init__(self, num_points=None, shape=None, points=None,
                 trim_domain=False, **kwargs):
        # Clean-up input points
        points = self._parse_points(shape=shape,
                                    num_points=num_points,
                                    points=points)
        # Create voronoi tessellation
        vor = sptl.Voronoi(points=points)
        # Convert to adjacency matrix
        am = topotools.vor_to_am(vor)
        # Extract conns
        conns = sp.vstack((am.row, am.col)).T
        # Make points back into 3D, if necessary
        points = vor.vertices
        if points.shape[1] == 2:
            points = sp.vstack((points.T, sp.zeros((points.shape[0], )))).T
        # Initialize network object
        super().__init__(Np=sp.shape(points)[0],
                         Nt=sp.shape(conns)[0],
                         **kwargs)
        # Update network with coords and conns from above
        self['pore.coords'] = sp.around(points, decimals=10)
        self['throat.conns'] = conns
        if trim_domain:
            Ps = self._find_external_pores(shape=shape)
            topotools.trim(network=self, pores=Ps)
