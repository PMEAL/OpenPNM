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

class GenericGeometry(OpenPNM.Base.Utilities):
    r"""
    GenericGeometry - Base class to construct pore networks

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
    network : OpenPNM Network Object
    
    name : string
        The name to apply to the labels (e.g. 'layer_1')
        
    locations : boolean mask or list of indices
        The pore locations in the network where this geometry applies
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages
        
    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> loc = pn.get_pore_indices() #Get all pores to define geometry everywhere
    >>> geo = OpenPNM.Geometry.GenericGeometry(name='geo_test',locations=loc,network=pn)
    >>> geo.add_method(prop='pore_seed',model='constant',value=0.123)
    >>> seeds = pn.get_pore_data(labels='geo_test',prop='pore_seed')
    >>> seeds[0]
    0.123
    """

    def __init__(self, network,name,locations,**kwargs):
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        loc = sp.array(locations,ndmin=1)
        if locations.dtype == bool:
            network.set_pore_info(label=name,locations=loc)
        else:
            network.set_pore_info(label=name,locations=loc)
        ind = network.get_pore_indices(name)
        r'''
        TODO: The following lines will create conflicting throat labels when additionaly geometries are added
        '''
        Tn = network.find_neighbor_throats(ind)
        network.set_throat_info(label=name,locations=Tn)
        network._geometry.append(self) #attach geometry to network
        self.name = name
        self._net = network #Attach network to self
        self._prop_list = []
              
    def regenerate(self):
        r'''
        This updates all properties using the selected methods
        '''
        self._logger.info("Refreshing geometry")
        for item in self._prop_list:
            self._logger.debug('Refreshing: '+item)
            getattr(self,item)()
    
    def add_method(self,prop='',**kwargs):
        try:
            function = getattr( getattr(OpenPNM.Geometry, prop), kwargs['model'] ) # this gets the method from the file
            preloaded_fn = partial(function, geometry=self, network=self._net, **kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])

if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    loc = pn.get_pore_indices()
    test = OpenPNM.Geometry.GenericGeometry(name='doc_test',locations=loc,network=pn)

