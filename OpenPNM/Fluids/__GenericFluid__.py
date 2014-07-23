"""
module Physics
===============================================================================

"""
import sys, os, collections
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
from functools import partial

class GenericFluid(OpenPNM.Utilities.Tools):
    r"""
    Base class to generate a generic fluid object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common fluids are included with OpenPNM and can be found under OpenPNM.Fluids.

    Parameters
    ----------
    network : OpenPNM Network object 
        The network to which this fluid should be attached
    name : str
        A unique string name to identify the fluid object, typically same as 
        instance name but can be anything.
    init_cond : dictionary, optional
        A dictionary of 'key':value pairs to initize fluid properties.  If omitted
        the temperature and pressure of the fluid are set to STP.  A typical 
        additiona would be 'mole_fraction' for a mixture.  Temperature and
        pressure can be set through init_cond, but will be set to STP if not 
        only other properties are included.
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    """
    def __init__(self,network,name=None,dynamic_data=False,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        # Attach objects to self for internal access
        self._net = network
        
        # Link this Fluid to the Network
        network._fluids.append(self) 
        
        # Initialize attributes
        self._physics = []
        self._models = collections.OrderedDict()
        self._dynamic_data = dynamic_data
        self.name = name
        
        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']        
        
        # Set default T and P since most propery models require it
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0
        
    def __getitem__(self,key):
        if key not in self.keys():
            return self.interleave_data(key)
        else:
            return super().__getitem__(key)
        
    def interleave_data(self,key):
        element = key.split('.')[0]
        temp = sp.ndarray((self.count(element),))
        for item in self._physics:
            locations = item.locations(element)
            if key not in item.keys():
                values = sp.ones_like(locations)*sp.nan
            else:
                values = item[key]
            temp[locations] = values
        return temp

    def regenerate(self,props=''):
        r'''
        This updates all properties of the fluid using the selected models
        
        Parameters
        ----------
        prop_list : string or list of strings
            The names of the properties that should be updated, defaults to all
        mode : string
            Control how the regeneration occurs.  
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> air = OpenPNM.Fluids.Air(loglevel=50,network=pn)
        >>> air.regenerate()  # Regenerate all properties at once
        >>> air.regenerate('molar_density')  # only one property
        >>> air.regenerate(['molar_density', 'diffusivity'])  # or several
        '''
        #First regenerate self
        if props == '':
            props = self._models.keys()
        elif type(props) == str:
            props = [props]
        for item in props:
            if item in self._models.keys():
                self[item] = self._models[item]()
            else:
                self._logger.warning('Requested proptery is not a dynamic model: '+item)
        
        #Then regenerate all physics objects associated with fluid
        for item in self._physics:
            item.regenerate()
            
        #Then pull in data from freshly regenerated Physics objects
        for phys in self._physics:
            for item in phys.props():
                element = item.split('.')[0]
                locations = self.locations(element)
                if item not in self.props():
                    self[item] = sp.nan
                self[item][locations] = phys[item]
        
    def add_model(self,model,propname,static=False,**kwargs):
        r'''
        Add specified property estimation model to the Fluid object.
        
        Parameters
        ----------
        na
        
        Examples
        --------
        None yet

        '''
        #Build partial function from given kwargs
        Ps = self.pores()
        Ts = self.throats()
        fn = partial(model,fluid=self,propname=propname,pores=Ps,throats=Ts,**kwargs)
        self[propname] = fn()  # Write static values to self
        if not static:  # Store model in a private ditionary
            self._models[propname] = fn
        
    def physics(self,name=''):
        r'''
        Retrieves Physics assocaiated with the Fluid
        
        Parameters
        ----------
        name : string, optional
            The name of the Physics object to retrieve
        Returns
        -------
            If name is NOT provided, then a list of Physics names is returned. 
            If a name IS provided, then the Physics object of that name is 
            returned.
        '''
        if name == '':
            phys = []
            for item in self._physics:
                phys.append(item.name)
        else:
            phys = self.find_object(obj_name=name)
        return phys

if __name__ =="__main__":
    pn = OpenPNM.Network.TestNet()
    fluid = OpenPNM.Fluids.GenericFluid(name='test_fluid',network=pn)
    import doctest
    doctest.testmod(verbose=True)



