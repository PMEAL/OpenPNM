# -*- coding: utf-8 -*-
"""
===============================================================================
module __Physics__: Base class for mananging pore-scale Physics properties
===============================================================================

"""
from OpenPNM.Base import logging, Tools
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

    def __init__(self, network=None, phase=None, geometry=None, pores=[],
                 throats=[], **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        # Associate with Network
        if network is None:
            network = GenericNetwork()
        self.network.update({network.name: network})  # Attach network to self
        self._net.physics.update({self.name: self})  # Register self with network

        # Associate with Phase
        if phase is None:
            phase = GenericPhase(network=self._net)
        phase.physics.update({self.name: self})  # Register self with phase
        self.phases.update({phase.name: phase})  # Register phase with self

        if geometry is not None:
            if (sp.size(pores) > 0) or (sp.size(throats) > 0):
                raise Exception('Cannot specify a Geometry AND pores or throats')
            pores = self._net.toindices(self._net['pore.' + geometry.name])
            throats = self._net.toindices(self._net['throat.' + geometry.name])

        # Initialize a label dictionary in the associated phase and network
        self._phases[0]['pore.'+self.name] = False
        self._phases[0]['throat.'+self.name] = False
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        try:
            self.set_locations(pores=pores, throats=throats)
        except:
            self.controller.purge_object(self)
            raise Exception('Provided locations are in use, instantiation cancelled')

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Convert self.name into 'all'
        if key.split('.')[-1] == self.name:
            key = element + '.all'
        if key in self.keys():  # Look for data on self...
            return super(GenericPhysics, self).__getitem__(key)
        else:  # ...Then check Network
            return self._phases[0][key][self._phases[0][element + '.' + self.name]]

    def _set_phase(self, phase):
        current_phase = self._phases[0]
        # Remove labels of self from current phase
        pore_label = current_phase.pop('pore.'+self.name)
        throat_label = current_phase.pop('throat.'+self.name)
        # Add labels of self to new phase
        phase['pore.'+self.name] = pore_label
        phase['throat.'+self.name] = throat_label
        # Replace phase reference on self
        self.phases.clear()
        self.phases.update({phase.name: phase})
        # Remove physics reference on current phase
        current_phase.physics.pop(self.name)
        phase.physics.update({self.name: self})

    def _get_phase(self):
        return self._phases[0]

    parent_phase = property(fget=_get_phase, fset=_set_phase)

    def set_locations(self, pores=[], throats=[], mode='add'):
        r"""
        Assign or unassign a Physics object to specified locations

        Parameters
        ----------
        pores : array_like
            The pore locations in the Network where this Physics is to apply

        throats : array_like
            The throat locations in the Network where this Physics is to apply

        mode : string
            Either 'add' (default) or 'remove' the object from the specified
            locations

        Examples
        --------
        >>> import OpenPNM

        """
        if mode == 'add':
            # Check if any constant values exist on the object
            for item in self.props():
                if (item not in self.models.keys()) or \
                   (self.models[item]['regen_mode'] == 'constant'):
                    raise Exception('Constant properties found on object, ' +
                                    'cannot increase size')
            if sp.size(pores) > 0:
                Tools.SetLocations.add(obj=self, element='pore',
                                       locations=pores)
            if sp.size(throats) > 0:
                Tools.SetLocations.add(obj=self, element='throat',
                                       locations=throats)
        if mode == 'remove':
            if sp.size(pores) > 0:
                Tools.SetLocations.drop(obj=self, element='pore',
                                        locations=pores)
            if throats is not None:
                Tools.SetLocations.drop(obj=self, element='throat',
                                        locations=throats)
        # Finally, regenerate models to correct the length of all arrays
        self.models.regenerate()
