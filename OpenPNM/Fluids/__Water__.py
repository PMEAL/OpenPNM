import OpenPNM
import scipy as sp

from __GenericFluid__ import GenericFluid

class Water(GenericFluid):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(Water,self).__init__(**kwargs)
        self._logger.debug("Construct class")

    def assign(self,network):
        params =       {'name': 'water',
                        'Pc': 2.206e6, #Pa
                        'Tc': 647,     #K
                        'MW': 0.0181,  #kg/mol
                        'diffusivity': {'method': 'constant',
                                        'DAB': 2e-12},
                        'viscosity': {'method': 'constant',
                                      'mu': 0.001},
                        'molar_density': {'method': 'constant',
                                          'c': 44444},
                        }
        network.phases.update({'water': params})

if __name__ =="__main__":
    testfluid_dict ={   'name':  'testfluid',
                       'phase':  1,
                 'diffusivity':  2.09e-5,
                   'viscosity':  1.73e-5,
                'molardensity':  40.9}
    pn = OpenPNM.Geometry.Cubic().generate()
    testfluid = GenericFluid(testfluid_dict,loggername="TestGenericFluid")
    testfluid.update(pn)