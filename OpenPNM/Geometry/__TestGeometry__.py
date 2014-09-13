"""
module __Toray090__: Subclass of GenericGeometry for a standard Toray TGPH090
gas diffusion layer.s
=============================================================================== 

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as _sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry.__GenericGeometry__ import GenericGeometry

class TestGeometry(GenericGeometry):
    r"""
    Toray090 subclass of GenericGeometry

    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super(TestGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
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
    pn = OpenPNM.Network.TestNet()
    pass
