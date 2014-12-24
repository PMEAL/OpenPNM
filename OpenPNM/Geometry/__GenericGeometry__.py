# -*- coding: utf-8 -*-
"""
===============================================================================
GenericGeometry -- Base class to manage pore scale geometry
===============================================================================

"""

import scipy as sp
from OpenPNM.Base import Core
from OpenPNM.Base import logging
from OpenPNM.Network import GenericNetwork
logger = logging.getLogger()
import OpenPNM.Geometry.models

class GenericGeometry(Core):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object

    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the physics nowhere.  The locations can be
        change after instantiation using ``set_locations()``.

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this this geometry applies.

    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> Ps = pn.pores()  # Get all pores
    >>> Ts = pn.throats()  # Get all throats
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts)
    """

    def __init__(self,network=None,pores=[],throats=[],seed=None,**kwargs):
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        logger.name = self.name

        if network is None:
            self._net = GenericNetwork()
        else:
            self._net = network  # Attach network to self
            self._net._geometries.append(self)  # Register self with network.geometries

        #Initialize a label dictionary in the associated network
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        self.set_locations(pores=pores,throats=throats)
        self._seed = seed

    def __getitem__(self,key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        else:
            return super(GenericGeometry,self).__getitem__(key)

    def set_locations(self,pores=[],throats=[]):
        r'''
        This method can be used to set the pore and throats locations of an
        *empty* object.  Once locations have been set they can not be changed.

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and/or throats where the object should be applied.

        Notes
        -----
        This method is intended to assist in the process of loading saved
        objects.  Save data can be loaded onto an empty object, then the object
        can be reassociated with a Network manually by setting the pore and
        throat locations on the object.
        '''
        pores = sp.array(pores,ndmin=1)
        throats = sp.array(throats,ndmin=1)
        if len(pores)>0:
            #Check for existing Geometry in pores
            temp = sp.zeros((self._net.Np,),bool)
            for key in self._net.geometries():
                temp += self._net['pore.'+key]
            overlaps = sp.sum(temp*self._net.tomask(pores=pores))
            if overlaps > 0:
                raise Exception('The given pores overlap with an existing Geometry object')
            #Initialize locations
            self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
            #Specify Geometry locations in Network dictionary
            self._net['pore.'+self.name][pores] = True
        if len(throats)>0:
            #Check for existing Geometry in pores
            temp = sp.zeros((self._net.Nt,),bool)
            for key in self._net.geometries():
                temp += self._net['throat.'+key]
            overlaps = sp.sum(temp*self._net.tomask(throats=throats))
            if overlaps > 0:
                raise Exception('The given throats overlap with an existing Geometry object')
            #Initialize locations
            self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
            #Specify Geometry locations in Network dictionary
            self._net['throat.'+self.name][throats] = True

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)

