
"""
module Physics
===============================================================================

"""
import OpenPNM
from functools import partial
import scipy as sp

class GenericPhysics(OpenPNM.Base.Utilities):
    
    def __init(self,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
    def create(self,network,fluid,**prms):
        r"""
        Create a fluid object using the supplied parameters
        """
        self.name = prms['name']
        self._prop_list = {}
        for key, args in prms.items():
            try:
                function = getattr( getattr(OpenPNM.Physics, key), args['method'] ) # this gets the method from the file
                preloaded_fn = partial(function, physics=self, network=network, fluid=fluid, **args) #
                setattr(self, key, preloaded_fn)
                self._logger.info("Successfully loaded {}.".format(key))
                self._prop_list[key] = True
            except AttributeError:
                self._logger.debug("Did not manage to load {}.".format(key))
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


