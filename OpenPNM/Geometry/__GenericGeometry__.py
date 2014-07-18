"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
from functools import partial

class GenericGeometry(OpenPNM.Utilities.Tools):
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
        self._models = {}
        
        #Initialize geometry locations
        self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
        self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
        network['pore.'+self.name] = False
        network['pore.'+self.name][pores] = True
        network['throat.'+self.name] = False
        network['throat.'+self.name][throats] = True
        
    def pores(self,**kwargs):
        return self._net.pores(labels=self.name)

    def throats(self,**kwargs):
        return self._net.throats(labels=self.name)

    def regenerate(self, props=''):
        r'''
        This updates all properties of the Geometry object
        
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
            prop_list = self.props()
        elif type(prop_list) == str:
            props = [prop_list]
        for item in prop_list:
            self[item] = self._models[item]()
        
    def add_model(self,model,propname,**kwargs):
        r'''
        Add specified pore scale model to the Geometry object
        
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
        if propname not in self._net.props():
            self._net[propname] = sp.nan
        self._net[propname][locations] = fn()
        self._models[propname] = fn
        self[propname] = fn()
        
    def geometry_health(self):
        geoms = self._net.geometries()
        temp = sp.zeros((self._net.Np,))
        for item in geoms:
            ind = self._net['pore.'+item]
            temp[ind] = temp[ind] + 1
        health = {}
        health['overlaps'] = sp.where(temp>1)[0].tolist()
        health['undefined'] = sp.where(temp==0)[0].tolist()
        return health

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)

