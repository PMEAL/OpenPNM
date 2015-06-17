# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericPhase__: Base class for building Phase objects
===============================================================================

"""
from OpenPNM.Network import GenericNetwork
import OpenPNM.Phases.models
from OpenPNM.Base import Core, Tools, logging
import scipy as sp
logger = logging.getLogger(__name__)


class GenericPhase(Core):
    r"""
    Base class to generate a generic phase object.  The user must specify models
    and parameters for all the properties they require. Classes for several
    common phases are included with OpenPNM and can be found under OpenPNM.Phases.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Phase should be attached

    components : list of OpenPNM Phase objects
        These Phase objects are ficticious or virtual phases that are the pure
        components from which the mixture is made.  They are used to calculate
        and store any pure component data.  If none are supplied then this
        object will act like either a pure component, a mixture whose properties
        are well known (like air) and need not to be found from consideration of
        the pure component properties.

    name : str, optional
        A unique string name to identify the Phase object, typically same as
        instance name but can be anything.

    """
    def __init__(self, network=None, components=[], **kwargs):
        super().__init__(**kwargs)
        logger.name = self.name

        if network is None:
            self._net = GenericNetwork()
        else:
            self._net = network

        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

        # Register Ohase object in Network dictionary
        self._net['pore.'+self.name] = True
        self._net['throat.'+self.name] = True

        if components != []:
            for comp in components:
                self.set_component(phase=comp)
        self._net._phases.append(self)  # Append this Phase to the Network

    def __setitem__(self, prop, value):
        for phys in self._physics:
            if (prop in phys.keys()) and ('all' not in prop.split('.')):
                    logger.error(prop + ' is already defined in at least one \
                                 associated Physics object')
                    return
        super().__setitem__(prop, value)

    def __getitem__(self, key):
        if key.split('.')[-1] == self.name:
            element = key.split('.')[0]
            return self[element+'.all']
        if key not in self.keys():
            logger.debug(key+' not on Phase, constructing data from Physics')
            return self._interleave_data(key, sources=self._physics)
        else:
            return super().__getitem__(key)

    def set_component(self, phase, mode='add'):
        r"""
        This method is used to add or remove a ficticious phase object to this
        object.

        Parameters
        ----------
        phase : OpenPNM Phase object
            This is the ficticious phase object defining a pure component.

        mode : string
            Indicates whether to 'add' or 'remove' the supplied Phase object
        """
        if mode == 'add':
            if phase.name in self.phases():
                logger.error('Phase already present')
                pass
            else:
                self._phases.append(phase)  # Associate any sub-phases with self
                phase._phases.append(self)  # Associate self with sub-phases
                # Add models for components to inherit mixture T and P
                phase.models.add(propname='pore.temperature',
                                 model=OpenPNM.Phases.models.misc.mixture_value)
                phase.models.add(propname='pore.pressure',
                                 model=OpenPNM.Phases.models.misc.mixture_value)
                # Move T and P models to beginning of regeneration order
                phase.models.reorder({'pore.temperature': 0, 'pore.pressure': 1})
        elif mode == 'remove':
            if phase.name in self.phases():
                self._phases.remove(phase)
                phase._phases = []
            else:
                logger.error('Phase not found')
                pass

    def check_mixture_health(self):
        r"""
        Query the properties of the 'virtual phases' that make up a mixture
        to ensure they all add up
        """
        mole_sum = sp.zeros((self.Np,))
        for comp in self._phases:
            try:
                mole_sum = mole_sum + comp['pore.mole_fraction']
            except:
                pass
        return mole_sum

    def check_physics_health(self):
        r"""
        Perform a check to find pores which have overlapping or undefined Physics
        """
        phys = self.physics()
        Ptemp = sp.zeros((self.Np,))
        Ttemp = sp.zeros((self.Nt,))
        for item in phys:
                Pind = self['pore.'+item]
                Tind = self['throat.'+item]
                Ptemp[Pind] = Ptemp[Pind] + 1
                Ttemp[Tind] = Ttemp[Tind] + 1
        health = Tools.HealthDict()
        health['overlapping_pores'] = sp.where(Ptemp > 1)[0].tolist()
        health['undefined_pores'] = sp.where(Ptemp == 0)[0].tolist()
        health['overlapping_throats'] = sp.where(Ttemp > 1)[0].tolist()
        health['undefined_throats'] = sp.where(Ttemp == 0)[0].tolist()
        return health
