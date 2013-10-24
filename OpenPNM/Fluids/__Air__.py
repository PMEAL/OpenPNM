import OpenPNM
import scipy as sp
import copy

from __GenericFluid__ import GenericFluid

class Air(GenericFluid):
    r"""
    Creates Fluid object with a default name 'air' and preset values
    """
    def __init__(self,**kwargs):
        super(Air,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
    def create(self,fluid_name='air'):
        r"""
        Creates Fluid object with a default name 'air'
        
        Parameters
        ----------

        fluid_name : string
            fluid_name = 'air' (default)\n
            fluid name that gets tagged to fluid-specific pore and throat conditions\n
        """
        self._fluid_recipe = {   'name': fluid_name,
                                   'Pc': 3.771e6, #Pa
                                   'Tc': 132.65,  #K
                                   'MW': 0.0291,  #kg/mol
                          'diffusivity': {'method': 'Fuller',
                                              'MA': 31.99,
                                              'MB': 28.01,
                                              'vA': 16.6,
                                              'vB': 17.9},
                            'viscosity': {'method': 'Reynolds',
                                              'uo': 0.001,
                                               'b': 0.1},
                        'molar_density': {'method': 'ideal_gas',
                                               'R': 8.314},
            }
        return self

if __name__ =="__main__":
    print 'test script not written for Air class'