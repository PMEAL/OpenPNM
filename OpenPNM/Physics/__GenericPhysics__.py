# -*- coding: utf-8 -*-
"""
===============================================================================
module __Physics__: Base class for mananging pore-scale Physics properties
===============================================================================

"""
from OpenPNM.Base import logging
from OpenPNM.Network import GenericNetwork
from OpenPNM.Phases import GenericPhase
import OpenPNM.Physics.models
import scipy as sp
logger = logging.getLogger(__name__)


class GenericPhysics(OpenPNM.Base.Core):
    r"""
    Generic class to generate Physics objects

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    geometry : OpenPNM Geometry object
        The Geometry object that defines the pores/throats where this Physics
        should be applied.  If this argument is supplied, then pores and
        throats cannot be specified.

    pores and/or throats : array_like
        The list of pores and throats where this physics applies. If either are
        left blank this will apply the physics nowhere.  The locations can be
        change after instantiation using ``set_locations()``.  If pores and
        throats are supplied, than a geometry cannot be specified.

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.

    """

    def __init__(self, network=None, phase=None, geometry=None,
                 pores=[], throats=[], **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        # Associate with Network
        if network is None:
            self._net = GenericNetwork()
        else:
            self._net = network  # Attach network to self
            self._net._physics.append(self)  # Register self with network

        # Associate with Phase
        if phase is None:
            self._phases.append(GenericPhase())
        else:
            phase._physics.append(self)  # Register self with phase
            self._phases.append(phase)  # Register phase with self

        if geometry is not None:
            if (sp.size(pores) > 0) or (sp.size(throats) > 0):
                raise Exception('Cannot specify a Geometry AND pores or throats')
            pores = self._net.toindices(self._net['pore.' + geometry.name])
            throats = self._net.toindices(self._net['throat.' + geometry.name])

        # Initialize a label dictionary in the associated fluid
        self._phases[0]['pore.'+self.name] = False
        self._phases[0]['throat.'+self.name] = False
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        self.set_locations(pores=pores, throats=throats)

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Convert self.name into 'all'
        if key.split('.')[-1] == self.name:
            key = element + '.all'

        if key in self.keys():  # Look for data on self...
            return super(GenericPhysics, self).__getitem__(key)
        else:  # ...Then check Network
            return self._phases[0][key][self._phases[0][element + '.' + self.name]]

    def set_locations(self, pores=[], throats=[], mode='add'):
        r"""
        Set the pore and throat locations of the Physics object

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and/or throats where the object should be applied.
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
