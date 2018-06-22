import scipy as sp
from openpnm.network import Cubic
from openpnm import topotools


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
