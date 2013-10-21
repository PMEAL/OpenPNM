import OpenPNM
import scipy as sp

class GenericFluid(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")

    def assign(self,network,params):
        network.phases.update({params['name']: params})

    def create(self,network,**params):
        self.fluid_dict = {}
        self.fluid_dict.update({'diffusivity': self.diffusivity(network,params['diffusivity'])})
        self.fluid_dict.update({'viscosity': self.viscosity(params['viscosity'])})
        self.fluid_dict.update({'molar_density': self.molar_density(params['molar_density'])})
        return self.fluid_dict

    def diffusivity(self,network,params):
        eqn = getattr(OpenPNM.Fluids.Diffusivity,params['method'])
        params.update({'network':network})
        vals = eqn(**params)
        del params['network']
        return vals

    def viscosity(self,params):
        eqn = getattr(OpenPNM.Fluids.Viscosity,params['method'])
        params['T'] = 50
        params['P'] = 20
        return eqn(**params)

    def molar_density(self,params):
        eqn = getattr(OpenPNM.Fluids.MolarDensity,params['method'])
        params['T'] = 50
        params['P'] = 20
        params['R'] = 8.314
        return eqn(**params)

if __name__ =="__main__":
    testfluid_dict ={   'name':  'testfluid',
                       'phase':  1,
                 'diffusivity':  2.09e-5,
                   'viscosity':  1.73e-5,
                'molardensity':  40.9}
    pn = OpenPNM.Geometry.Cubic().generate()
    testfluid = GenericFluid(testfluid_dict,loggername="TestGenericFluid")
    testfluid.update(pn)