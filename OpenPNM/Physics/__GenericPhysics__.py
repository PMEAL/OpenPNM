
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
        
    def create(self,net,fluid,**prms):
        r"""
        Create a fluid object using the supplied parameters
        """
        self.name = prms['name']
        for key, args in prms.items():
            try:
                function = getattr( getattr(OpenPNM.Physics, key), args['method'] ) # this gets the method from the file
                preloaded_fn = partial(function, physics=self, network=net, fluid=fluid, **args) #
                setattr(self, key, preloaded_fn)
                self._logger.info("Successfully loaded {}.".format(key))
            except AttributeError:
                self._logger.debug("Did not manage to load {}.".format(key))
        self.regenerate()
        return self

    def regenerate(self):
        r'''
        This updates all properties using the methods indicated in the recipe.  This method also takes the opportunity to ensure all values are Numpy arrays.
        '''
        try: self.capillary_pressure()
        except: pass
        try: self.thermal_conductance()
        except: pass
        try: self.hydraulic_conductance()
        except: pass
        try: self.diffusive_conductance()
        except: pass
        try: self.electronic_conductance()
        except: pass


