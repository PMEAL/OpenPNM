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
        self._implemented_methods = ['diffusivity','viscosity','molar_density']
        self._accepted_throat_methods = []

    def create(self,fluid_dict):
        self._fluid_dict = copy.deepcopy(fluid_dict)
        return self
        

    def assign_to_network(self,network):
        #Make sure base conditions are set
        if 'temperature' not in network.pore_conditions.keys():
            network.pore_conditions['temperature'] = sp.array(293.15,ndmin=1)
        if 'pressure' not in network.pore_conditions.keys():
            network.pore_conditions['pressure'] = sp.array(101325,ndmin=1)
        network.phases.update({self._fluid_dict['name']: self._fluid_dict})
        self.refresh_in_network(network)

    def remove_from_network(self,network):
        fluid_name = self._fluid_dict['name']
        del network.phases[fluid_name]
        for i in self._implemented_methods:
            try: del network.pore_conditions["{prop}_{fluid}".format(prop=i, fluid=fluid_name)]
            except: self._logger.warn("Attempted to remove a pore condition that doesn't exist")
#        for i in self._accepted_throat_methods:
#            try: del network.throat_conditions["{prop}_{fluid}".format(prop=i, fluid=fluid_name)]
#            except: self._logger.warn("Attempted to remove a throat condition that doesn't exist")

    def rename_fluid(self,new_name):
        self._fluid_dict['name'] = new_name
        self._logger.warn("Warning: renaming a fluid does not change fluids already assigned to a network")

    def refresh_in_network(self,network):
        fluid_name = self._fluid_dict['name']
        try: network.phases[fluid_name]
        except: raise Exception('This fluid does not exist in the specified network')
        for i in self._implemented_methods:
            network.pore_conditions.update({ "{prop}_{fluid}".format(prop=i, fluid=fluid_name) : getattr(self,i)(network)})
#        for i in self._accepted_throat_methods:
#            network.throat_conditions.update({ "{prop}_{fluid}".format(prop=i, fluid=fluid_name) : getattr(self,i)(network)})

    def diffusivity(self,network):
        fluid_name = self._fluid_dict['name']
        params = network.phases[fluid_name]['diffusivity']
        eqn = getattr(OpenPNM.Fluids.Diffusivity,params['method'])
        vals = eqn(network,**params)
        return sp.array(vals,ndmin=1)

    def viscosity(self,network):
        fluid_name = self._fluid_dict['name']
        params = network.phases[fluid_name]['viscosity']
        eqn = getattr(OpenPNM.Fluids.Viscosity,params['method'])
        vals = eqn(network,**params)
        return sp.array(vals,ndmin=1)

    def molar_density(self,network):
        fluid_name = self._fluid_dict['name']
        params = network.phases[fluid_name]['molar_density']
        eqn = getattr(OpenPNM.Fluids.MolarDensity,params['method'])
        vals = eqn(network,**params)
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
