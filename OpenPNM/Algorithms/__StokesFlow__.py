"""
===============================================================================
module __StokesFlow__: Viscous fluid flow
===============================================================================

"""
import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class StokesFlow(GenericLinearTransport):
    r'''
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the network.
        
    '''
    
    def __init__(self,**kwargs):
        r'''
        '''
        super(StokesFlow,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,conductance='hydraulic_conductance',quantity='pressure',**params):
        r'''
        '''
        self._logger.info("Setup "+self.__class__.__name__)         
        super(StokesFlow,self).setup(conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()
                
    def calc_eff_permeability(self):
        r'''
        '''
        D_normal = self._calc_eff_prop()
        self._eff_property = D_normal*sp.mean(self._phase['pore.viscosity'])
        return self._eff_property
        