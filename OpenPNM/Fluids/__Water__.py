
from .__GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r"""
    Creates Fluid object with a default name 'water' and preset values
    """
    def __init__(self,**kwargs):
        super(Water,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.recipe = {
        'Name': 'water',
        'Pc': 2.206e6, #Pa
        'Tc': 647,     #K
        'MW': 0.0181,  #kg/mol
        'diffusivity': {'method': 'constant',
                        'value': 1e-12},
        'viscosity': {'method': 'constant',
                      'value': 0.001},
        'molar_density': {'method': 'constant',
                          'value': 44445},
        'surface_tension': {'method': 'Eotvos',
                            'k': 2.25e-4},
        'contact_angle': {'method': 'constant',
                          'value': 120},
        }

if __name__ =="__main__":
    print('no tests yet')
