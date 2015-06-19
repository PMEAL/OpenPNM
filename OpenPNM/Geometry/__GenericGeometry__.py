# -*- coding: utf-8 -*-
"""
===============================================================================
GenericGeometry -- Base class to manage pore scale geometry
===============================================================================

"""

import scipy as sp
import OpenPNM.Geometry.models
from OpenPNM.Base import Core
from OpenPNM.Base import logging
from OpenPNM.Network import GenericNetwork
logger = logging.getLogger(__name__)


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
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    """

    def __init__(self, network=None, pores=[], throats=[], seed=None, **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        if network is None:
            self._net = GenericNetwork()
        else:
            self._net = network  # Attach network to self
            # Register self with network.geometries
            self._net._geometries.append(self)

        # Initialize a label dictionary in the associated network
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        self.set_locations(pores=pores, throats=throats)
        self._seed = seed

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Convert self.name into 'all'
        if key.split('.')[-1] == self.name:
            key = element + '.all'

        if key in self.keys():  # Look for data on self...
            return super(GenericGeometry, self).__getitem__(key)
        if key == 'throat.conns':  # Handle specifically
            [P1, P2] = \
                self._net['throat.conns'][self._net[element + '.' + self.name]].T
            Pmap = sp.zeros((self._net.Np,), dtype=int) - 1
            Pmap[self._net.pores(self.name)] = self.Ps
            conns = sp.array([Pmap[P1], Pmap[P2]]).T
            # Replace -1's with nans
            if sp.any(conns == -1):
                conns = sp.array(conns, dtype=object)
                conns[sp.where(conns == -1)] = sp.nan
            return conns
        else:  # ...Then check Network
            return self._net[key][self._net[element+'.'+self.name]]

    def set_locations(self, pores=[], throats=[], mode='add'):
        r"""
        Set the pore and throat locations of the Geometry object

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and/or throats in the Network where the object
            should be applied
        mode : string
            Indicates whether list of pores or throats is to be added or removed
            from the object.  Options are 'add' (default) or 'remove'.

        """
        if len(pores) > 0:
            pores = sp.array(pores, ndmin=1)
            self._set_locations(element='pore', locations=pores, mode=mode)
        if len(throats) > 0:
            throats = sp.array(throats, ndmin=1)
            self._set_locations(element='throat', locations=throats, mode=mode)
