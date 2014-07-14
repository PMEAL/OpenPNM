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
    
    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the physics nowhere.  The locations can be
        change after instantiation using ``set_locations()``.
    
    name : str, optional
        A unique string name to identify the Physics object, typically same as 
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.  
    
    """

    def __init__(self,network,fluid,pores=[],throats=[],name=None,**kwargs):
        super(GenericPhysics,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        
        #Setup containers for ojecct linking
        self._prop_list = []
        
        # Append objects for internal access
        self._net = network
        self._fluid = fluid
        
        self.name = name
        
        #Use composition to assign pores and throats to this physics
        fluid['pore.'+self.name] = False
        fluid['throat.'+self.name] = False
        fluid['pore.'+self.name][pores] = True
        fluid['throat.'+self.name][throats] = True
        fluid['pore.'+self.name] = False
        fluid['throat.'+self.name] = False
        fluid['pore.'+self.name][pores] = True
        fluid['throat.'+self.name][throats] = True
        self.pores = partial(fluid.pores,labels=self.name)
        self.throats = partial(fluid.throats,labels=self.name)
        
    @property
    def Np(self):
        return self.num_pores()
        
    @property
    def Nt(self):
        return self.num_throats()
    
    def num_pores(self):
        return sp.shape(self.pores())[0]
        
    def num_throats(self):
        return sp.shape(self.throats())[0]
        
    def count(self,element):
        #Remove pluralizaton
        if element == 'pores':
            element = 'pore'
        if element == 'throats':
            element = 'throat'
        temp = {}
        temp['pore'] = self.num_pores()
        temp['throat'] = self.num_throats()
        if element != None:
            temp = temp[element]
        return temp
        
    def generate(self):
        raise NotImplementedError('This method must be implemented in a subclass')
    
    def set_locations(self,pores=[],throats=[],mode='add'):
        r'''
        Assign Physics object to specifed pores and/or throats
        '''
        if pores != []:
            if mode == 'add':
                self._fluid['pore.'+self.name][pores] = True
            elif mode == 'overwrite':
                self._fluid['pore.'+self.name] = False
                self._fluid['pore.'+self.name][pores] = True
            elif mode == 'remove':
                self._fluid['pore.'+self.name][pores] = False
            else:
                print('invalid mode received')
        if throats != []:
            if mode == 'add':
                self._fluid['throat.'+self.name][throats] = True
            elif mode == 'overwrite':
                self._fluid['throat.'+self.name] = False
                self._fluid['throat.'+self.name][throats] = True
            elif mode == 'remove':
                self._fluid['throat.'+self.name][throats] = False
            else:
                print('invalid mode received')
        
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
            prop_attr = item.replace('.','_')
            element = 'throat'
            prop_key = element+'.'+prop_attr.split(element+'_')[-1]
            if element == 'pore':
                locations = self.pores()
            elif element == 'throat':
                locations = self.throats()
            if prop_key not in self._fluid.keys():
                self._fluid[prop_key] = sp.nan
            self._logger.debug('Refreshing: '+item)
            self._fluid[prop_key][locations] = getattr(self,prop_attr)()
            
    def add_method(self,prop='',prop_name='',**kwargs):
        r'''
        THIS METHOD IS DEPRECATED USE add_property() INSTEAD
        '''
        self.add_property(prop=prop,prop_name=prop_name,**kwargs)
            
    def add_model(self,model,propname,**kwargs):
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
        if propname not in self._fluid.keys():
            self._fluid[propname] = sp.ones((self.count(element),))*sp.nan
        self._fluid[propname][fn.keywords[locations]] = fn()
        
    def add_property(self,prop='',prop_name='',**kwargs):
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
            function = getattr( getattr(OpenPNM.Physics.models, prop), kwargs['model'] ) # this gets the method from the file
            if prop_name: prop = prop_name #overwrite the default prop with user supplied name  
            preloaded_fn = partial(function, network=self._net, fluid=self._fluid, pores=self.pores(),throats=self.throats(),**kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])

if __name__ == '__main__':
    print('none yet')


