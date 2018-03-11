# -*- coding: utf-8 -*-
"""
===============================================================================
CubicTemplate: Generate cubic lattice networks of arbitrary overall shape
===============================================================================

"""
import scipy as sp
from openpnm.network import Cubic
from openpnm import topotools

shapes = [[10, 10, 10], [20, 20, 5], [5, 10, 8]]
spacings = [1, 0.5, [2, 1, 0.8]]
origins = [[0, 0, 0], [0, 0, 10], [0, 0, 12.5]]
labels = None


class CubicTemplate(Cubic):
    r"""
    """
    def __init__(self, template, **kwargs):

        template = sp.atleast_3d(template)
        super().__init__(shape=template.shape, **kwargs)

        coords = sp.unravel_index(range(template.size), template.shape)
        self['pore.template_coords'] = sp.vstack(coords).T
        self['pore.template_indices'] = self.Ps
        self['pore.drop'] = template.flatten() == 0
        topotools.trim(network=self, pores=self.pores('drop'))
        del self['pore.drop']

    def toarray(self, prop):
        pass

    def fromarray(self, array):
        pass
