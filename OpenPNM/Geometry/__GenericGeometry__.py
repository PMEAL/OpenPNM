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

class GenericGeometry(OpenPNM.Utilities.Base,dict):
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
    >>> Ps = pn.pores()  # Get all pores
    >>> Ts = pn.throats()  # Get all throats
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn)
    >>> geom.set_locations(pores=Ps,throats=Ts)
    """

    def __init__(self,network,pores=[],throats=[],name=None,**kwargs):
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        self._net = network #Attach network to self        
        self.name = name
        #Register self with network.geometries
        self._net._geometries.append(self)
        
        #Initialize geometry locations if given
        network['pore.'+self.name] = False
        network['pore.'+self.name][pores] = True
        network['throat.'+self.name] = False
        network['throat.'+self.name][throats] = True
        self.num_pores = partial(network.num_pores,labels=self.name)
        self.num_throats = partial(network.num_throats,labels=self.name)
        
    def generate(self):
        raise NotImplementedError('This method must be implemented in a subclass')
        
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
            if mode == 'add':
                self._net.set_info(label=self.name,pores=pores,mode='merge')
            elif mode == 'remove':
                self._net.set_info(label=self.name,pores=pores,mode='remove')
            else:
                print('invalid mode received')
        if throats != []:
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
    
    def regenerate(self, props=''):
        r'''
        This updates all properties using the selected methods
        
        Parameters
        ----------
        prop_list : string or list of strings
            The names of the properties that should be updated, defaults to all
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> Ps = pn.pores()
        >>> geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,pores=Ps)
        >>> geom.generate()
        >>> geom.regenerate()  # Regenerate all properties at once
        >>> geom.regenerate('pore.seed')  # only one property
        >>> geom.regenerate(['pore.seed', 'pore.diameter'])  # or several
        '''
        if props == '':
            prop_list = self.keys()
        elif type(prop_list) == str:
            props = [prop_list]
        for item in prop_list:
            element = item.split('.')[0]
            locations = self[item].keywords[element+'s']
            self._net[item][locations] = self[item]()
        
    def add_model(self,model,propname,**kwargs):
        r'''
        Add specified pore scale model to the Geometry object.
        
        Parameters
        ----------
        model : function
            The model retrieved from the ./Geometry/models library
        propname : string
            The name of the physical property calculated by the model.  This 
            name is used as the dictionary key in the Network object.
        kwargs : keyword arguments
            These are the arguments required by the model, with the exception
            of network, pores and throats which are passed automatically.
        
        Examples
        --------
        None yet

        '''
        #Build partial function from given kwargs
        fn = partial(model,network=self._net,propname=propname,pores=self.pores(),throats=self.throats(),**kwargs)
        #Determine element and locations
        element = propname.split('.')[0]
        locations = fn.keywords[element+'s']
        if propname not in self._net.keys():
            self._net[propname] = sp.nan
        self._net[propname][locations] = fn()
        self[propname] = fn

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)

