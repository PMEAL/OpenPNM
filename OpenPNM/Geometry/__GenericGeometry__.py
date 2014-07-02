"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM
import scipy as sp
from functools import partial

class GenericGeometry(OpenPNM.Utilities.Base):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object
    
    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this this geometry applies.
    
    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)

    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages
        
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> loc = pn.pores() #Get all pores to define geometry everywhere
    >>> geo = OpenPNM.Geometry.GenericGeometry(network=pn)
    >>> geo.set_locations(pores=loc)
    >>> geo.add_method(prop='pore_seed',model='constant',value=0.123)
    >>> geo.regenerate()
    >>> seeds = pn.get_pore_data(locations=geo.name,prop='seed')
    >>> seeds[0]
    0.123
    """

    def __init__(self,network,name=None,**kwargs):
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        self._net = network #Attach network to self        
        self.name = name
        #Setup containers for object linking
        self._physics = [] #Create list for physics to append themselves to
        self._net._geometries.append(self)
        self._prop_list = []
        
        #Initialize geometry to NOWHERE
        network.set_info(label=self.name,pores=[])
        network.set_info(label=self.name,throats=[])
        self.num_pores = partial(network.num_pores,labels=self.name)
        self.num_throats = partial(network.num_throats,labels=self.name)
        
    @property
    def Np(self):
        return self.num_pores()
        
    @property
    def Nt(self):
        return self.num_throats()
        
    def count(self,element=None):
        temp = {}
        temp['pore'] = self.num_pores()
        temp['throat'] = self.num_throats()
        if element != None:
            temp = temp[element]
        return temp
    
    def set_locations(self,pores=[],throats=[],mode='add'):
        r'''
        Assign Geometry object to specifed pores (or throats)
        '''
        if pores != []:
            if pores == 'all':
                pores = self._net.pores(labels='all')
            if mode == 'add':
                self._net.set_info(label=self.name,pores=pores,mode='merge')
            elif mode == 'remove':
                self._net.set_info(label=self.name,pores=pores,mode='remove')
            else:
                print('invalid mode received')
        
        if throats != []:
            if throats == 'all':
                throats = self._net.get_throat_indices(labels='all')
            if mode == 'add':
                self._net.set_info(label=self.name,throats=throats,mode='merge')
            elif mode == 'remove':
                self._net.set_info(label=self.name,throats=throats,mode='remove')
            else:
                print('invalid mode received')
        
    def pores(self):
        r'''
        '''
        return self._net.pores(labels=self.name)

    def throats(self):
        r'''
        '''
        return self._net.throats(labels=self.name)
    
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
        >>> pn = OpenPNM.Network.TestNet()
        >>> pind = pn.pores()
        >>> geom = OpenPNM.Geometry.Stick_and_Ball(network=pn)
        >>> geom.regenerate()  # Regenerate all properties at once
        >>> geom.regenerate('pore_seed')  # only one property
        >>> geom.regenerate(['pore_seed', 'pore_diameter'])  # or several
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
        None yet
        '''
        try:
            function = getattr( getattr(OpenPNM.Geometry, prop), kwargs['model'] ) # this gets the method from the file
            if prop_name: propname = prop = prop_name #overwrite the default prop with user supplied name
            else:
                #remove leading pore_ or throat_ from dictionary key
                propname = prop.split('_')[1]
                element = prop.split('_')[0]
                if len(prop.split('_')) > 2:
                    propname = prop.split(element+'_')[1] 
            preloaded_fn = partial(function, geometry=self, network=self._net,propname=propname, **kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])

    def check_consistency(self):
        r'''
        Checks to see if the current geometry conflicts with other geometries, 
        and ensures that geometries have been applied to all locations
        '''
        temp = sp.zeros_like(self._net.get_pore_info(label=self.name),dtype=int)
        for item in self._net._geometries:
            temp = temp + sp.array(self._net.get_pore_info(label=item),dtype=int)
        print('Geometry labels overlap in', sp.sum(temp>1),'pores')
        print('Geometry not yet applied to',sp.sum(temp==0),'pores')
        
        temp = sp.zeros_like(self._net.get_throat_info(label=self.name),dtype=int)
        for item in self._net._geometries:
            temp = temp + sp.array(self._net.get_throat_info(label=item),dtype=int)
        print('Geometry labels overlap in', sp.sum(temp>1),'throats')
        print('Geometry not yet applied to',sp.sum(temp==0),'throats')

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)

