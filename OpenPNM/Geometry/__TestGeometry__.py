# -*- coding: utf-8 -*-
"""
===============================================================================
TestGeometry -- A Geomery for Toray TGPH090 gas diffusion layers
===============================================================================

"""

import scipy as _sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry

class TestGeometry(GenericGeometry):
    r"""
    Toray090 subclass of GenericGeometry

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(TestGeometry,self).__init__(**kwargs)
        self._generate()

    def _generate(self):
        r'''
        '''
        self.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       seed=1)
        self.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        self['pore.diameter'] = self['pore.seed']
        self['throat.diameter'] = self['throat.seed']
        self['pore.volume']=_sp.pi/6*self['pore.diameter']**3
        self['pore.area']=_sp.constants.pi/4*self['pore.diameter']**2

        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self['throat.volume'] = _sp.pi/4*self['throat.length']*self['throat.diameter']**2
        self['throat.area'] = _sp.constants.pi/4*(self['throat.diameter'])**2
        self['throat.surface_area'] = _sp.constants.pi/(self['throat.diameter'])*self['throat.length']

if __name__ == '__main__':
    import OpenPNM
    pn = OpenPNM.Network.TestNet()
    pass
