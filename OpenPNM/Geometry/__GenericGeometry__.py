"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
import scipy.stats as spst
from functools import partial
import numpy as np

class GenericGeometry(OpenPNM.Base.Utilities):
    r"""
    GenericGeometry - Base class to construct pore networks

    This class contains the interface definition for the construction of networks

    Parameters
    ----------
    network : OpenPNM Network Object
    
    name : string
        The name to apply to the subdomain (e.g. 'layer_1')
        
    locations : boolean mask or list of indices
        The pore locations in the network where this geometry applies
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    """

    def __init__(self, network,name,locations,**kwargs):
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Method: Constructor")
        loc = sp.array(locations,ndmin=1)
        if locations.dtype == bool:
            network.set_pore_info(prop=name,locations=loc,is_indices=False)
        else:
            network.set_pore_info(prop=name,locations=loc,is_indices=True)
        ind = network.get_pore_indices(name)
        r'''
        TODO: The following lines will create conflicting throat labels when additionaly geometries are added
        '''
        Tn = network.get_neighbor_throats(ind)
        network.set_throat_info(prop=name,locations=Tn,is_indices=True)
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
    pn = OpenPNM.Network.Cubic(name='test_net').generate(divisions=[5,5,5],lattice_spacing=[0.001])
    test = OpenPNM.Geometry.GenericGeometry(loglevel=10,name='doc_test',locations=[0],network=pn)

