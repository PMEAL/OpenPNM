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
    def __init__(self,network,name=None,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        # Attach objects for internal access
        self._net = network
        
        # Link this Fluid to the Network
        network._fluids.append(self) 
        
        # Initialize tracking lists
        self._physics = []
        self._prop_list = []
        
        self.name = name
        
        # Set default T and P since most propery models require it
        self.set_data(prop='temperature',pores='all',data=298.0)
        self.set_data(prop='pressure',pores='all',data=101325.0)
        
        # Initialize label 'all' in the object's own info dictionaries
        self.set_info(label='all',pores=network.pores('all'))
        self.set_info(label='all',throats=network.throats('all'))
        
    def apply_conditions(self,**values):
        r'''
        Apply multiple scalar conditions to the fluid in a single step
        
        Parameters
        ----------
        values : name / value pair ()
            Any arguments can be passed, along with a corresponding values.  
            These arguments will be assigned to the fluid's pore data
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> air = OpenPNM.Fluids.Air(loglevel=50,network=pn)
        >>> air.apply_conditions(molar_mass=18,density=1000)
        '''
        for item in values.keys():
            self.set_data(prop=item,pores='all',data=values[item])

    def regenerate(self,prop_list='',mode=None):
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
        
    def add_property(self,prop='',prop_name='',**kwargs):
        r'''
        Add specified property estimation model to the fluid object.
        
        Parameters
        ----------
        prop : string
            The name of the fluid property attribute to add.
            This name must correspond with a file in the Fluids folder.  
            To add a new property simply add a file with the appropriate name and the necessary methods.
        prop_name : string, optional
            This argument will be used as the method name and the dictionary key
            where data is written by method. This option is provided for occasions
            when multiple properties of the same type are required, such as
            diffusivity coefficients of each species in a multicomponent mixture.
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> fluid = OpenPNM.Fluids.GenericFluid(network=pn)
        >>> fluid.add_method(prop='diffusivity',model='constant',value=1.234)
        >>> fluid.regenerate()
        >>> fluid.get_pore_data(prop='diffusivity')
        array([ 1.234])
        
        In cases where the same property is needed multiple times (like 
        diffusivity of multiple species in a mixture), it is possible to
        specify unique property names:
        
        >>> fluid.add_method(prop='diffusivity',prop_name='diffusivity_of_species_2',model='constant',value=3.333)
        >>> fluid.regenerate()
        '''
        try:
            function = getattr( getattr(OpenPNM.Fluids, prop), kwargs['model'] ) # this gets the method from the file
            if prop_name: prop = prop_name #overwrite the default prop with user supplied name  
            preloaded_fn = partial(function, fluid=self, network=self._net, propname=prop, **kwargs) #          
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])
        
    def physics_listing(self):
        r"""
        Prints the names of all physics objects attached to the fluid
        """
        for item in self._physics:
            print(item.name+': ',item)

    def physics_update(self,name='all'):
        r"""
        Updates ALL properties of specified physics object attached to the fluid

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
    pn = OpenPNM.Network.TestNet()
    fluid = OpenPNM.Fluids.GenericFluid(name='test_fluid',network=pn)
    import doctest
    doctest.testmod(verbose=True)



