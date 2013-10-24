import OpenPNM
import scipy as sp
import copy

from __GenericFluid__ import GenericFluid

class Toray090(GenericFluid):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(Toray090,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        # the following two lines track the methods that update pore or throat properties
        self._implemented_methods = ['diffusivity','viscosity','molar_density']
        self._accepted_throat_methods = []

    def create(self,**ignored):
        self._fluid_recipe =   {'name': 'Toray090',
                                        'thermal_conductivity': {'method': 'constant',
                                                                 'DAB': 2e-12},
                               }
        return self

    def ThermalConductivity(self,fluid):
        r"""
        In-plane thermal conductivity of Carbon paper without Teflon based on Zamel et al.  
    
        """
        T = fluid.pore_conditions['temperature'] -273.15    
        k = -7.166e-6*T**3 + 2.24e-3*T**2-0.237*T + 20.1
        return sp.array(k,ndmin=1)

if __name__ =="__main__":
    testfluid_dict ={   'name':  'testfluid',
                       'phase':  1,
                 'diffusivity':  2.09e-5,
                   'viscosity':  1.73e-5,
                'molardensity':  40.9}
    pn = OpenPNM.Geometry.Cubic().generate()
    testfluid = GenericFluid(testfluid_dict,loggername="TestGenericFluid")
    testfluid.update(pn)