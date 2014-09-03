"""
===============================================================================
module __OhmicConduction__: Electronic or ionic conduction
===============================================================================

"""
import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class OhmicConduction(GenericLinearTransport):
    r'''
    A subclass of GenericLinearTransport to simulate electron and ionic 
    conduction.  The 2 main roles of this subclass are to set the default 
    property names and to implement a method for calculating the effective 
    conductivity of the network.

    '''

    def __init__(self,**kwargs):
        r'''
        '''
        super(OhmicConduction,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,conductance='electronic_conductance',quantity='voltage',**params):
        r'''
        '''  
        self._logger.info("Setup "+self.__class__.__name__)        
        super(OhmicConduction,self).setup(conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()