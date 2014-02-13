"""
module Physics
===============================================================================

"""
import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)

import OpenPNM
import scipy as sp
from functools import partial

class GenericPhysics(OpenPNM.Base.Utilities):
    r"""
    Base class to generate a generic Physics object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common Physics are included with OpenPNM and can be found under OpenPNM.Physics.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this Physics should be attached
        
    fluid : OpenPNM Fluid object 
        The Fluid object to which the Physics applies
    
    name : str
        A unique string name to identify the Physics object, typically same as 
        instance name but can be anything.
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    """
    def __init__(self,network,fluid,name,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.name = name
        self._prop_list = []
        self._fluid = []
        #bind objects togoether
        self._fluid.append(fluid) #attach fluid to physics
        fluid._physics.append(self)
        self._net = network

    def regenerate(self):
        r'''
        This updates all properties using the selected methods
        '''
        for item in self._prop_list:
            self._logger.debug('Refreshing: '+item)
            getattr(self,item)()
            
    def add_method(self,prop='',prop_name='',**kwargs):
        r'''
        Add specified property estimation model to the physics object.
        
        Parameters
        ----------
        prop : string
            The name of the pore scale physics property attribute to add.
            This name must correspond with a file in the Physics folder.  
            To add a new property simply add a file with the appropriate name and the necessary methods.
           
        prop_name : string, optional
            This argument will be used as the method name and the dictionary key
            where data is written by method. This option is provided for occasions
            when multiple properties of the same type are required, such as
            diffusive conductance of each species in a multicomponent mixture.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        '''
        try:
            function = getattr( getattr(OpenPNM.Physics, prop), kwargs['model'] ) # this gets the method from the file
            if prop_name: prop = prop_name #overwrite the default prop with user supplied name  
            preloaded_fn = partial(function, physics=self, network=self._net, propname=prop, fluid=self._fluid[0], **kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])

if __name__ == '__main__':
    print('none yet')


