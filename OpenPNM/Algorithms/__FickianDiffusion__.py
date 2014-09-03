"""
===============================================================================
module __FickianDiffusion__: Diffusive mass transfer
===============================================================================

"""
import scipy as sp
from .__GenericLinearTransport__ import GenericLinearTransport

class FickianDiffusion(GenericLinearTransport):
    r'''
    A subclass of GenericLinearTransport to simulate binary diffusion.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the effective diffusion coefficient
    of the network.
    '''

    def __init__(self,**kwargs):
        r'''
        
        '''
        super(FickianDiffusion,self).__init__(**kwargs)
        self._logger.info('Create '+self.__class__.__name__+' Object')
        
    def run(self,conductance='diffusive_conductance',quantity='mole_fraction',**params):
        r'''
        '''  
        self._logger.info("Setup "+self.__class__.__name__)   
        super(FickianDiffusion,self).setup(conductance=conductance,quantity=quantity)
        
        super(GenericLinearTransport,self).run()
        
    def calc_eff_diffusivity(self):
        r'''
        '''
        D_normal = self._calc_eff_prop()
        self._eff_property = D_normal/sp.mean(self._phase['pore.molar_density'])
        return self._eff_property
        



