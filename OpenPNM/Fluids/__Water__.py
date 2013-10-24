import OpenPNM
import scipy as sp
import copy

from __GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r"""
    Creates Fluid object named 'water'
    """
    def __init__(self,**kwargs):
        super(Water,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
    def create(self,name='water'):
        self._fluid_recipe = {   'name': name,
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