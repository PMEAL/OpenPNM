"""

module __Permeability__: 
===============================================================================

"""

import scipy as sp
from .__LinearSolver__ import LinearSolver

class StokesFlow(LinearSolver):
    r'''
        
    '''
    
    def __init__(self,**kwargs):
        r'''
        '''
        super(StokesFlow,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def setup(self,fluid,conductance='hydraulic_conductance',quantity='pressure',**params):
        r'''
        '''
        self._logger.info("Setup "+self.__class__.__name__)         
        super(StokesFlow,self).setup(fluid=fluid,conductance=conductance,quantity=quantity)
                
    def calc_eff_permeability(self, clean=False):
        D_normal = self._calc_eff_prop()
        self._eff_property = D_normal*self._fluid['pore.viscosity']
        return self._eff_property
        