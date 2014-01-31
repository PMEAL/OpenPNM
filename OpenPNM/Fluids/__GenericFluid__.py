
"""
module Physics
===============================================================================

"""
import OpenPNM
import scipy as sp
from functools import partial

class GenericFluid(OpenPNM.Base.Tools):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.pore_data = {}
        self.throat_data = {}
        self.pore_info = {}
        self.throat_info = {}
        self._physics = []

    def create(self,network,T=298.,P=101325.,**recipe):
        r"""
        Create a fluid object using the supplied parameters
        """
        try: recipe = self.recipe #check if recipe is pre-existing on self (from init of subclassed methods)
        except: pass
        try: self.name = recipe['name']
        except: self._logger.error('Fluid name must be given')
        #Bind objects together
        network._fluids.append(self)
        self.set_pore_data(prop='numbering',data=network.get_pore_indices()) #This is necessary for the methods from 'tools' to work.  They must know network size.
        self.set_throat_data(prop='numbering',data=network.get_throat_indices())        
        self.Tc = recipe['Tc']
        self.Pc = recipe['Pc']
        self.MW = recipe['MW']
        self.set_pore_data(prop='temperature',data=sp.array(T,ndmin=1))
        self.set_pore_data(prop='pressure',data=sp.array(P,ndmin=1))
        for key, args in recipe.items():
            try:
                function = getattr( getattr(OpenPNM.Fluids, key), args['method'] ) #Get method from the file
                preloaded_fn = partial(function, fluid=self, network=network, **args) 
                setattr(self, key, preloaded_fn)
                self._logger.info('Successfully added '+key+' to '+self.name)
            except AttributeError: pass
        self.regenerate()
        return self

    def regenerate(self):
        r'''
        This updates all properties using the methods indicated in the recipe.
        '''
        try: self.viscosity()
        except: pass
        try: self.diffusivity()
        except: pass
        try: self.molar_density()
        except: pass
        try: self.surface_tension()
        except: pass
        try: self.contact_angle()
        except: pass
        #Update physics associated with this fluid too
        self.physics_update()
        
    def physics_listing(self):
        r"""
        Prints the names of all physics objects attached to the network
        """
        for item in self._physics:
            print(item.name+': ',item)

    def physics_update(self,name='all'):
        r"""
        Updates ALL properties of specified physics object attached to the network

        Parameters
        ----------
        name : string (optional)
            The name of physics object to be updated.  An empty string (default) refreshes all physics.
        """
        for item in self._physics:
            if (item.name == name) or (name == 'all'):
                item.regenerate()
                self._logger.info('Refreshed '+item.name)
        
    def __str__(self):
        return('This is the __str__ methods of the generic_fluid being overwritten')

if __name__ =="__main__":

    #Create fluids
    air_recipe = {
    'name': 'air',
    'Pc': 3.771e6, #Pa
    'Tc': 132.65,  #K
    'MW': 0.0291,  #kg/mol
    'diffusivity': {'method': 'Fuller',
                    'MA': 0.03199,
                    'MB': 0.0291,
                    'vA': 16.3,
                    'vB': 19.7},
    'viscosity': {'method': 'Reynolds',
                  'uo': 0.001,
                  'b': 0.1},
    'molar_density': {'method': 'ideal_gas',
                      'R': 8.314},
    'surface_tension': {'method': 'constant',
                        'value': 0},
    'contact_angle': {'method': 'na'},
    }
    gas = OpenPNM.Fluids.GenericGas().create(air_recipe)



