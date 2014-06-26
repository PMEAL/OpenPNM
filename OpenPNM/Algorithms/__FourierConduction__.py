"""

module __FourierConduction__
===============================================================================

"""

import scipy as sp
from .__LinearSolver__ import LinearSolver

class FourierConduction(LinearSolver):
    r"""
    
    """
    
    def __init__(self,**kwargs):
        r'''
        '''
        super(FourierConduction,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def setup(self,fluid,conductance='thermal_conductance',quantity='temperature',**params):
        r'''
        '''  
        self._logger.info('Setup '+self.__class__.__name__)         
        super(FourierConduction,self).setup(fluid=fluid,conductance=conductance,quantity=quantity)
        