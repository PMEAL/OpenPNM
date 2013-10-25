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
        #List of fluid property categories that are invoked when fluid is created
        self._implemented_methods = ['diffusivity',
                                     'viscosity',
                                     'molar_density',
                                     'surface_tension',
                                     'contact_angle']

    def create(self,fluid_recipe):
        self._fluid_recipe = copy.deepcopy(fluid_recipe)
        self.pore_conditions = {}
        self.throat_conditions = {}
        self.pore_conditions.update({'temperature': 298})
        self.pore_conditions.update({'pressure': 101325})
        self.regenerate()
        return self

    def refresh(self):
        self.regenerate()

    def regenerate(self):
        r'''
        This updates all properties using the methods indicated in the recipe.  This method also takes the opportunity to ensure all values are Numpy arrays.
        '''
        for condition in self._implemented_methods:
            self.pore_conditions.update({condition: getattr(self,condition)()})
        #Make sure all values are Numpy arrays (Np,) or (Nt,)
        for i in self.pore_conditions.keys():
            self.pore_conditions[i] = sp.array(self.pore_conditions[i],ndmin=1)
        for i in self.throat_conditions.keys():
            self.throat_conditions[i] = sp.array(self.throat_conditions[i],ndmin=1)

    def reset(self):
        r'''
        Remove all existing condtions from the fluid

        TODO: This works, but is kludgy
        '''
        self.pore_conditions = {}
        self.throat_conditions = {}
        try: del self.partner
        except: pass
        self.pore_conditions.update({'temperature': 298})
        self.pore_conditions.update({'pressure': 101325})
        self.regenerate()
        return self

    def clone(self):
        r'''
        Create an exact duplicate fluid, but a unique object.

        TODO: Doesn't work yet
        '''
        return self.__new__

    def set_pair(self,fluid2):
        r'''
        Creates a fluid pair by storing each fluid object on the other.  This allows tracking of defending vs invading phase, among other things.

        TODO: This method needs plenty of checks, for preexisting pair, etc.

        '''
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

    def contact_angle(self):
        params = self._fluid_recipe['contact_angle']
        eqn = getattr(OpenPNM.Fluids.ContactAngle,params['method'])
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
