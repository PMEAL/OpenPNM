
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
    def __init__(self,network,name,T=298,P=101325,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.name = name
        self.pore_data = {}
        self.throat_data = {}
        self.pore_info = {}
        self.throat_info = {}
        self._physics = []
        self._prop_list = []
        self._net = network
        self.set_pore_data(prop='temperature',data=T)
        self.set_pore_data(prop='pressure',data=P)

    def regenerate(self):
        r'''
        This updates all properties using the selected methods
        '''
        for item in self._prop_list:
            self._logger.debug('Refreshing: '+item)
            getattr(self,item)()
        
    def add_method(self,prop='',**kwargs):
        try:
            function = getattr( getattr(OpenPNM.Fluids, prop), kwargs['model'] ) # this gets the method from the file
            preloaded_fn = partial(function, fluid=self, network=self._net, **kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])
        
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



