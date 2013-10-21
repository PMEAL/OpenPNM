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

    def create(self,name='water',**params):
        self.fluid_dict = {}
        self.fluid_dict['name'] = name
        self.fluid_dict.update({'diffusivity': self.diffusivity()})
        self.fluid_dict.update({'viscosity': self.viscosity()})
        self.fluid_dict.update({'molar_density': self.molar_density()})
        return self.fluid_dict

    def update(self,**params):
        self.create(**params)

    def diffusivity(self,**params):
        return 0.0

    def viscosity(self,**params):
        return 0.001

    def molar_density(self,**params):
        return 1000

if __name__ =="__main__":
    testfluid_dict ={   'name':  'testfluid',
                       'phase':  1,
                 'diffusivity':  2.09e-5,
                   'viscosity':  1.73e-5,
                'molardensity':  40.9}
    pn = OpenPNM.Geometry.Cubic().generate()
    testfluid = GenericFluid(testfluid_dict,loggername="TestGenericFluid")
    testfluid.update(pn)