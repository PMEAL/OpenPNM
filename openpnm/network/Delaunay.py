"""
===============================================================================
Delaunay: Generate random networks based on the Delaunay Tessellation
===============================================================================

"""
import scipy as sp
import scipy.spatial as sptl
from openpnm import topotools
from openpnm.network import GenericTessellation
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Delaunay(GenericTessellation):
    r"""

    """

    def __init__(self, num_points=None, shape=None, points=None,
                 trim_domain=True, **kwargs):
        # Clean-up input points
        points = self._parse_points(shape=shape,
                                    points=points,
                                    num_points=num_points)
        # Perform Delaunay tessellation on points
        tri = sptl.Delaunay(points)
        # Convert tessellation to adjacency matrix
        am = topotools.tri_to_am(tri)
        # Extract conns
        conns = sp.vstack((am.row, am.col)).T
        # Make points 3D if necessary
        if points.shape[1] == 2:
            points = sp.vstack((points.T, sp.zeros((points.shape[0], )))).T
        # Initialize network object now that Np and Nt are known
        super().__init__(Np=points.shape[0], Nt=conns.shape[0], **kwargs)
        # Assign conns and coords
        self['throat.conns'] = conns
        self['pore.coords'] = points
        if trim_domain:
            Ps = self._find_external_pores(shape=shape)
            topotools.trim(network=self, pores=Ps)
        # Attach delaunay object in case someone needs it?
        self._tri = tri
