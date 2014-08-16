"""

module __FourierConduction__
===============================================================================

"""

import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class FourierConduction(GenericLinearTransport):
    r"""
    
    """
    
    def __init__(self,**kwargs):
        r'''
        '''
        super(FourierConduction,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,phase,conductance='thermal_conductance',quantity='temperature',**params):
        r'''
        '''  
        self._logger.info('Setup '+self.__class__.__name__)         
        super(FourierConduction,self).setup(phase=phase,conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()