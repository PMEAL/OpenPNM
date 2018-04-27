"""
===============================================================================
Gabriel: Generate random networks based on the Gabriel Tessellation
===============================================================================

"""
from openpnm.network import Delaunay
from openpnm import topotools
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Gabriel(Delaunay):
    r"""

    """
    def __init__(self, num_points=None, shape=None, points=None, **kwargs):

        super().__init__(num_points=num_points, shape=shape, points=points,
                         **kwargs)
        topotools.trim(self, throats=self.throats('gabriel', mode='complement'))
        del self['throat.gabriel']
