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

class GenericPhysics(OpenPNM.Utilities.Base):
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

    def __init__(self,network,fluid,geometry,name='default',**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        self.name = name
        
        #Setup containers for ojecct linking
        self._prop_list = []

        # Append objects for internal access
        self._net = network
        self._fluid = fluid
        self._geometry = geometry
        
        # Connect this physics with it's geometry
        geometry._physics.update({name:self})
        fluid._physics.update({name:self})

        #Use composition to assign pores and throats to this physics
        self.pores = geometry.pores
        self.throats = geometry.throats

    def regenerate(self, prop_list='',mode=None):
        r'''
        This updates all properties using the selected methods

        Parameters
        ----------
        prop_list : string or list of strings
            The names of the properties that should be updated, defaults to all
        mode : string
            Control how the regeneration occurs.  
            
        Examples
        --------
        For examples refer to usage of Fluid or Geometry refresh methods
        '''
        if prop_list == '':
            prop_list = self._prop_list
        elif type(prop_list) == str:
            prop_list = [prop_list]
        if mode == 'exclude':
            a = sp.array(self._prop_list)
            b = sp.array(prop_list)
            c = a[sp.where(~sp.in1d(a,b))[0]]
            prop_list = list(c)
        for item in prop_list:
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
            preloaded_fn = partial(function, physics=self, network=self._net, propname=prop, fluid=self._fluid, geometry=self._geometry, **kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])

if __name__ == '__main__':
    print('none yet')


