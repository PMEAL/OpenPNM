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
from OpenPNM.Base import Tools

class GenericFluid(Tools):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,network,name,T=298,P=101325,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.name = name
        self._net = network
        #Initialize necessary empty attributes
        self.pore_data = {}
        self.throat_data = {}
        self.pore_info = {}
        self.throat_info = {}
        self._physics = []
        self._prop_list = []
        #Set default T and P since most propery models require it
        self.set_pore_data(prop='temperature',data=T)
        self.set_pore_data(prop='pressure',data=P)
        #Initialize 'numbering arrays in the objects own info dictionaries
        self.set_pore_info(label='all',locations=self._net.get_pore_indices())
        self.set_throat_info(label='all',locations=self._net.get_throat_indices())

    def regenerate(self):
        r'''
        This updates all properties using the selected methods
        '''
        for item in self._prop_list:
            self._logger.debug('Refreshing: '+item)
            getattr(self,item)()
        
    def add_method(self,prop='',**kwargs):
        r'''
        Add specified property estimation method to the fluid object.
        
        Parameters
        ----------
        prop : string
            The name of the fluid property attribute to add.
            This name must correspond with a file in the Fluids folder.  
            To add a new property simply add a file with the appropriate name and the necessary methods.
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> print(pn.name)
        test_network
        >>> fluid = OpenPNM.Fluids.GenericFluid(network=pn,name='test_fluid')
        >>> fluid.add_method(prop='diffusivity',model='constant',value=1.234)
        >>> fluid.regenerate()
        >>> fluid.get_pore_data(prop='diffusivity') #Use fluid's getter
        array([ 1.234])
        >>> pn.get_pore_data(prop='diffusivity',phase=fluid) #Use network's getter
        array([ 1.234])
        '''
        try:
            function = getattr( getattr(OpenPNM.Fluids, prop), kwargs['model'] ) # this gets the method from the file
            preloaded_fn = partial(function, fluid=self, network=self._net, **kwargs) #
            setattr(self, prop, preloaded_fn)
            self._logger.info("Successfully loaded {}.".format(prop))
            self._prop_list.append(prop)
        except AttributeError: print('could not find',kwargs['model'])
        
    def physics_listing(self):
        r"""
        Prints the names of all physics objects attached to the network
        """
        for item in self._physics:
            print(item.name+': ',item)

    def physics_update(self,name='all'):
        r"""
        Updates ALL properties of specified physics object attached to the network

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



