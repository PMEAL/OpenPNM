import OpenPNM
import scipy as sp
import copy

class GenericFluid(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        # the following two lines track the methods that update pore or throat properties
        self._implemented_methods = ['diffusivity','viscosity','molar_density','surface_tension']
        self._accepted_throat_methods = []
        self.pore_conditions = {}
        self.throat_conditions = {}

    def create(self,fluid_recipe):
        self._fluid_recipe = copy.deepcopy(fluid_recipe)
        return self

    def refresh(self):
        for condition in self._implemented_methods:
            self.pore_conditions.update({condition: getattr(self,condition)()})

    def regenerate(self):
        for condition in self._implemented_methods:
            self.pore_conditions.update({condition: getattr(self,condition)()})

    def set_pair(self,fluid2):
        self.partner = fluid2
        fluid2.partner = self

    def diffusivity(self):
        params = self._fluid_recipe['diffusivity']
        eqn = getattr(OpenPNM.Fluids.Diffusivity,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def viscosity(self):
        params = self._fluid_recipe['viscosity']
        eqn = getattr(OpenPNM.Fluids.Viscosity,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def molar_density(self):
        params = self._fluid_recipe['molar_density']
        eqn = getattr(OpenPNM.Fluids.MolarDensity,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def surface_tension(self):
        params = self._fluid_recipe['surface_tension']
        eqn = getattr(OpenPNM.Fluids.SurfaceTension,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

if __name__ =="__main__":

    pn = OpenPNM.Geometry.Cubic(loglevel=40).generate()

    #Define the fluids and set their properties
    params_air = {       'name': 'air',
                           'Pc': 3.771e6, #Pa
                           'Tc': 132.65,  #K
                           'MW': 0.0291,  #kg/mol
                  'diffusivity': {'method': 'Fuller',
                                      'MA': 31.99,
                                      'MB': 28.01,
                                      'vA': 16.6,
                                      'vB': 17.9},
                    'viscosity': {'method': 'Reynolds',
                                      'uo': 0.001,
                                       'b': 0.1},
                'molar_density': {'method': 'ideal_gas',
                                       'R': 8.413},
                   }
    #Create fluids
    air = OpenPNM.Fluids.GenericFluid(params_air)

    #Assign fluids to network
    air.assign(pn)
    print ''
    print 'current pore conditions:'
    for i in pn.pore_conditions.keys():
        print i,'=',pn.pore_conditions[i]
    print ''
