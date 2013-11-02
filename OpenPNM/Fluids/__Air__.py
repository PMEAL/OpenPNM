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

    def create(self,fluid_name='air',T=298.,P=101325.):
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
                                              'MA': 0.03199,
                                              'MB': 0.0291,
                                              'vA': 16.3,
                                              'vB': 19.7},
                            'viscosity': {'method': 'Reynolds',
                                              'uo': 0.001,
                                               'b': 0.1},
                        'molar_density': {'method': 'ideal_gas',
                                               'R': 8.314},
                      'surface_tension': {'method': 'na'},
                        'contact_angle': {'method': 'na'},
                        }
        self.pore_conditions = {}
        self.throat_conditions = {}
        self.pore_conditions.update({'temperature': T})
        self.pore_conditions.update({'pressure': P})
        self.regenerate()
        return self

if __name__ =="__main__":
    print 'test script not written for Air class'