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
        The Fluid object to which this Physics applies
        
    geometry : OpenPNM Geometry object
        The Geometry object to which this Physics applies
    
    name : str
        A unique string name to identify the Physics object, typically same as 
        instance name but can be anything.
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    """

    def __init__(self,network,fluid,geometry,name=None,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        #Setup containers for ojecct linking
        self._prop_list = []

        # Append objects for internal access
        self._net = network
        self._fluid = fluid
        self._geometry = geometry
        
        # Connect this physics with it's geometry
        geometry._physics.append(self)
        fluid._physics.append(self)

        self.name = name

        #Use composition to assign pores and throats to this physics
        self.pores = geometry.pores
        self.throats = geometry.throats
        self.Np = geometry.Np
        self.Nt = geometry.Nt
        self.count = geometry.count

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
        THIS METHOD IS DEPRECATED USE add_property() INSTEAD
        '''
        self.add_property(prop=prop,prop_name=prop_name,**kwargs)
            
    def add_property(self,model,propname,**kwargs):
        r'''
        Add specified property estimation model to the fluid object.
        
        Parameters
        ----------
        na
        
        Examples
        --------
        None yet

        '''
        #Determine element and locations
        element = propname.split('.')[0]
        if element == 'pore':
            locations = 'pores'
        elif element == 'throat':
            locations = 'throats'
        #Build partial function from given and updated kwargs
        self._fluid
        self._net
        self.pores()
        self.throats()
        fn = partial(model,fluid=self._fluid,network=self._net,pores=self.pores(),throats=self.throats(),**kwargs)
        if propname not in self._net.keys():
            self._fluid[propname] = sp.ones((self.count(element),))*sp.nan
        self._fluid[propname][fn.keywords[locations]] = fn()

if __name__ == '__main__':
    print('none yet')


