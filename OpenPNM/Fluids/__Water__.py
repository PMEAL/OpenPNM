import scipy as sp
import copy

from .__GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r"""
    Creates Fluid object with a default name 'water' and preset values
    """
    def __init__(self,**kwargs):
        super(Water,self).__init__(**kwargs)
        self._logger.debug("Construct class")

    def create(self,T=298.,P=101325.):
        r"""
        Creates Fluid object with a default name 'water'

        Parameters
        ----------

        fluid_name : string
            fluid_name = 'water' (default)\n
            fluid name that gets tagged to fluid-specific pore and throat conditions\n
        """
        water_recipe = {
        'Name': 'water',
        'Thermo':   {'Pc': 2.206e6, #Pa
                     'Tc': 647,     #K
                     'MW': 0.0181,  #kg/mol
                     },
        'Diffusivity': {'method': 'constant',
                        'value': 1e-12},
        'Viscosity': {'method': 'constant',
                      'value': 0.001},
        'MolarDensity': {'method': 'constant',
                          'value': 44445},
        'SurfaceTension': {'method': 'Eotvos',
                            'k': 2.25e-4},
        'ContactAngle': {'method': 'constant',
                          'value': 120},
        }
        self.pore_conditions = {}
        self.throat_conditions = {}
        self.pore_conditions.update({'temperature': T})
        self.pore_conditions.update({'pressure': P})
        self.regenerate()
        return self

if __name__ =="__main__":
    print('no tests yet')
