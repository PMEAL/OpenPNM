"""
===============================================================================
Gabriel: Generate random networks based on the Gabriel Tessellation
===============================================================================

"""
import scipy as sp
import scipy.spatial as sptl
from openpnm.network import Delaunay
from openpnm.topotools import trim
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Gabriel(Delaunay):
    r"""

    """
    def __init__(self, **kwargs):
        # Generate Delaunay tessellation from super class, then trim
        super().__init__(**kwargs)
        points = self['pore.coords']
        conns = self['throat.conns']
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
        trim(self, throats=~g)
