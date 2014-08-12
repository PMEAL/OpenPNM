"""

module __FickianDiffusion__
===============================================================================

"""

import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class FickianDiffusion(GenericLinearTransport):
    r'''
    '''

    def __init__(self,**kwargs):
        r'''
        '''
        super(FickianDiffusion,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,fluid,conductance='diffusive_conductance',quantity='mole_fraction',**params):
        r'''
        '''  
        self._logger.info("Setup "+self.__class__.__name__)   
        super(FickianDiffusion,self).setup(fluid=fluid,conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()
        
    def calc_eff_diffusivity(self, clean=False):
        D_normal = self._calc_eff_prop()
        self._eff_property = D_normal/sp.mean(self._fluid['pore.molar_density'])
        return self._eff_property
        



