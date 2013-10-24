import OpenPNM
import scipy as sp
import copy

from __GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r"""
    Creates Fluid object with a default name 'water' and preset values
    """
    def __init__(self,**kwargs):
        super(Water,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
    def create(self,fluid_name='water'):
        r"""
        Creates Fluid object with a default name 'water'
        
        Parameters
        ----------

        fluid_name : string
            fluid_name = 'water' (default)\n
            fluid name that gets tagged to fluid-specific pore and throat conditions\n
        """
        self._fluid_recipe = {   'name': fluid_name,
                                   'Pc': 2.206e6, #Pa
                                   'Tc': 647,     #K
                                   'MW': 0.0181,  #kg/mol
                          'diffusivity': {'method': 'constant',
                                           'value': 1e-12},
                            'viscosity': {'method': 'constant',
                                           'value': 0.001},
                        'molar_density': {'method': 'constant',
                                           'value': 44445},
                           }
        return self

if __name__ =="__main__":
    print 'test script not written for Water class'