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

    def create(self,params):
        self.fluid_dict = params
        return self.fluid_dict

    def assign(self,network,params):
        network.phases.update({params['name']: params})
        return self.refresh(network,params['name'])

    def remove(self,network,name):
        del network.phases[name]

    def refresh(self,network,name):
        params = network.phases[name]
        print params
        self.fluid_dict = {}
        self.fluid_dict.update({'diffusivity': self.diffusivity(network,params['diffusivity'])})
        self.fluid_dict.update({'viscosity': self.viscosity(network,params['viscosity'])})
        self.fluid_dict.update({'molar_density': self.molar_density(network,params['molar_density'])})
        for i in self.fluid_dict.keys():
            network.pore_conditions.update({i+'_'+name: self.fluid_dict[i]})
        return self.fluid_dict
        
    def diffusivity(self,network,params):
        eqn = getattr(OpenPNM.Fluids.Diffusivity,params['method'])
        vals = eqn(network,**params)
        return sp.array(vals,ndmin=1)

    def viscosity(self,network,params):
        eqn = getattr(OpenPNM.Fluids.Viscosity,params['method'])
        vals = eqn(network,**params)
        return sp.array(vals,ndmin=1)

    def molar_density(self,network,params):
        eqn = getattr(OpenPNM.Fluids.MolarDensity,params['method'])
        vals = eqn(network,**params)
        return sp.array(vals,ndmin=1)

if __name__ =="__main__":
    testfluid_dict ={   'name':  'testfluid',
                       'phase':  1,
                 'diffusivity':  2.09e-5,
                   'viscosity':  1.73e-5,
                'molardensity':  40.9}
    pn = OpenPNM.Geometry.Cubic().generate()
    testfluid = GenericFluid(testfluid_dict,loggername="TestGenericFluid")
    testfluid.update(pn)