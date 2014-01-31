
"""
module Physics
===============================================================================

"""
import OpenPNM
import scipy as sp
from functools import partial

class GenericPhysics(OpenPNM.Base.Utilities):
    r"""
    GenericPhysics - Base class to generate pore scale physics properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self._prop_list = {}
        self._fluid = []
        
    def create(self,network,fluid,**recipe):
        r"""
        Create a fluid object using the supplied parameters
        """
        try: recipe = self.recipe #check if recipe is pre-existing on self (from init of subclassed methods)
        except: pass
        try: self.name = recipe['name']
        except: self._logger.error('Physics name must be given')
        #bind objects togoether
        self._fluid.append(fluid) #attach fluid to physics
        fluid._physics.append(self)
        for key, args in recipe.items():
            try:
                function = getattr( getattr(OpenPNM.Physics, key), args['method'] ) # this gets the method from the file
                preloaded_fn = partial(function, physics=self, network=network, fluid=fluid, **args) #
                setattr(self, key, preloaded_fn)
                self._logger.info("Successfully loaded {}.".format(key))
                self._prop_list[key] = True
            except AttributeError: pass
        self.regenerate()
        return self

    def regenerate(self):
        r'''
        This updates all properties using the methods indicated in the recipe.
        '''
        for item in self._prop_list.keys():
            getattr(self,item)()
#        try: self.capillary_pressure()
#        except: self._logger.error('Failed to load')
#        try: self.thermal_conductance()
#        except: self._logger.error('Failed to load')
#        try: self.hydraulic_conductance()
#        except: self._logger.error('Failed to load')
#        try: self.diffusive_conductance()
#        except: self._logger.error('Failed to load')
#        try: self.electronic_conductance()
#        except: self._logger.error('Failed to load')


