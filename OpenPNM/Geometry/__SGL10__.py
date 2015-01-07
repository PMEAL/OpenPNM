# -*- coding: utf-8 -*-
"""
===============================================================================
SGL10 -- A geometry model for SGL10 type Gas Diffusion Layers
===============================================================================

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry

class SGL10(GenericGeometry):
    r"""
    SGL10 subclass of GenericGeometry.

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(SGL10,self).__init__(**kwargs)
        self._generate()

    def _generate(self):
        r'''
        '''
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       num_range=[0,0.8834],
                       seed=self._seed)
        self.remove_model('pore.seed')
        self.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        self.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=3.07,
                       psd_loc=1.97e-6,
                       psd_scale=1.6e-5,
                       psd_offset=18e-6)
        self.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        self.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        self.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=3.07,
                       tsd_loc=1.97e-6,
                       tsd_scale=1.6e-5,
                       tsd_offset=18e-6)
        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)


if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    pass