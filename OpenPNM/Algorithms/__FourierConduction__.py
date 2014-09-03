"""
===============================================================================
module __FourierConduction__: Conductive heat transfer
===============================================================================

A subclass of GenericLinearTransport to simulate heat conduction

"""

import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class FourierConduction(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate heat conduction.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the effective conductivity of the network.
    """
    
    def __init__(self,**kwargs):
        r'''
        '''
        super(FourierConduction,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,conductance='thermal_conductance',quantity='temperature',**params):
        r'''
        '''  
        self._logger.info('Setup '+self.__class__.__name__)         
        super(FourierConduction,self).setup(conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()