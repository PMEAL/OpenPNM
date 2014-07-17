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

class GenericPhysics(OpenPNM.Utilities.Base,dict):
    r"""
    Generic class to generate Physics objects  

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
        
        #Append objects for internal access
        self._net = network
        self._fluid = fluid
        self._net._physics.append(self)
        self._fluid._physics.append(self)
        self._models = {}
        
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
        
    def count(self,element=None):
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
        
    def regenerate(self, props=''):
        r'''
        This updates all properties using the selected methods

        Parameters
        ----------
        props : string or list of strings
            The names of the properties that should be updated, defaults to all
            
        Examples
        --------
        na
        '''
        if props == '':
            prop_list = self.keys()
        elif type(prop_list) == str:
            props = [prop_list]
        for item in prop_list:
            self[item] = self._models[item]()
            
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
        fn = partial(model,fluid=self._fluid,network=self._net,pores=self.pores(),throats=self.throats(),**kwargs)
        if propname not in self._fluid.keys():
            self._fluid[propname] = sp.ones((self.count(element),))*sp.nan
        self._fluid[propname][fn.keywords[locations]] = fn()
        self._models[propname] = fn
        self[propname] = fn()
        
    def fluids(self):
        temp = []
        temp.append(self._fluid.name)
        return temp
        
if __name__ == '__main__':
    print('none yet')


