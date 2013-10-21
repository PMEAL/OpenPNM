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

    def create(self,fluid):
        self.fluid_dict = fluid
        return self.fluid_dict

    def assign(self,network,fluid):
        network.phases.update({fluid['name']: fluid})
        self.refresh(network,fluid['name'])

    def remove(self,network,fluid_name):
        del network.phases[fluid_name]

    def rename(self,fluid,fluid_name):
        fluid['name'] = fluid_name

    def refresh(self,network,fluid_name):
        params = network.phases[fluid_name]
        self.fluid_dict = {}
        self.fluid_dict.update({'diffusivity': self.diffusivity(network,params['diffusivity'])})
        self.fluid_dict.update({'viscosity': self.viscosity(network,params['viscosity'])})
        self.fluid_dict.update({'molar_density': self.molar_density(network,params['molar_density'])})
        for i in self.fluid_dict.keys():
            network.pore_conditions.update({ "{prop}_{fluid}".format(prop=i, fluid=fluid_name) : self.fluid_dict[i]})
        
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