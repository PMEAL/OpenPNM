from openpnm.phases import GenericPhase as GenericPhase
import openpnm.models.phases as _models
import numpy as np
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class MultiPhase(GenericPhase):
    r"""
    Creates Phase object that represents a multiphase system consisting of
    a given list of OpenPNM Phase objects.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    project : OpenPNM Project object, optional
        The Project with which this phase should be associted.  If a
        ``network`` is given then this is ignored and the Network's project
        is used.  If a ``network`` is not given then this is mandatory.

    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    Examples
    --------
    >>> import scipy as sp
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> air = op.phases.Air(network=pn, name='air')  # Create two pure phases
    >>> water = op.phases.Water(network=pn, name='water')
    >>> m = op.phases.MultiPhase(network=pn, phases=[air, water], name='multi')
    >>> Ps = pn['pore.coords'][:, 0] < 3  # Pick some pores to be air filled
    >>> Ts = pn.find_neighbor_throats(pores=Ps)  # Find neighboring throats
    >>> Ts = pn.tomask(throats=Ts)  # Convert throat indices to mask
    >>> m.set_occupancy(phase=air, Pvals=Ps, Tvals=Ts)  # Assign occupancies
    >>> m.set_occupancy(phase=water, Pvals=~Ps, Tvals=~Ts)

    Air and water have uniform viscosity values throughout, so the peak-to-peak
    distance is 0, while the mixture phase has the viscosity of air in some
    locations and water in others, hence a heterogenous viscosity array:

    >>> sp.ptp(water['pore.viscosity']) == 0
    True
    >>> sp.ptp(air['pore.viscosity']) == 0
    True
    >>> sp.ptp(m['pore.viscosity']) == 0
    False

    """
    def __init__(self, phases=[], settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update({'phases': []})
        self.settings.update(settings)

        self['pore.occupancy.all'] = np.zeros(self.Np, dtype=float)
        self['throat.occupancy.all'] = np.zeros(self.Nt, dtype=float)

        self.pop('pore.temperature', None)
        self.pop('pore.pressure', None)

        # Add any supplied phases to the phases list
        for phase in phases:
            self.settings['phases'].append(phase.name)

        logger.warning('MultiPhases are a beta feature and functionality may' +
                       ' change in future versions')

    def __getitem__(self, key):
        try:
            vals = super().__getitem__(key)
        except KeyError:
            vals = self.interleave_data(key)
        return vals

    def _update_occupancy(self):
        # Update "occupancy.all" by summing up all occupancies
        for elem in ["pore", "throat"]:
            dict_ = self[f"{elem}.occupancy"]
            dict_.pop(f"{elem}.occupancy.all")
            self[f"{elem}.occupancy.all"] = np.sum(list(dict_.values()), axis=0)


    def _get_phases(self):
        phases = {self.project[item].name: self.project[item]
                  for item in self.settings['phases']}
        return phases

    phases = property(fget=_get_phases)

    def interleave_data(self, prop):
        r"""
        Gathers property values from component phases to build a single array

        If the requested ``prop`` is not on this MultiPhase, then a search is
        conducted on all associated phase objects, and values from each
        are assembled into a single array.

        Parameters
        ----------
        prop : string
            The property to be retrieved

        Returns
        -------
        array : ND-array
            An array containing the specified property retrieved from each
            component phase and assembled based on the specified mixing rule

        """
        element, _ = prop.split('.')
        vals = np.zeros([self._count(element=element)], dtype=float)
        # Retrieve property from constituent phases (weight = occupancy)
        try:
            for phase in self.phases.values():
                vals += phase[prop] * self[element + '.occupancy.' + phase.name]
        # Otherwise - if not found - retrieve from super class
        except KeyError:
            vals = super().interleave_data(prop)

        # Check for consistency of occypancy values (i.e. add up to 1)
        if np.any(self[element + '.occupancy.all'] != 1.0):
            self._update_occupancy()
            if np.any(self[element + '.occupancy.all'] != 1.0):
                raise Exception(f"Occupancy doesn't add to unity in all {element}s")
        return vals

    def regenerate_models(self, **kwargs):
        r"""
        Regenerate models associated with the Multiphase object

        This method works by first regenerating the models associated with the
        constituent phases, and then regenerating Multiphase models.

        """
        # Regenerate models associated with phases within MultiPhase object
        for phase in self.phases.values():
            phase.regenerate_models(**kwargs)
        # Regenerate models specific to MultiPhase object
        super().regenerate_models(**kwargs)

    def set_occupancy(self, phase, Pvals=[], Tvals=[], pores=[], throats=[]):
        r"""
        Specify occupancy of a phase in each pore and/or throat

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase whose occupancy is being specified

        Pvals : array_like, float
            The volume fraction of ``phase`` in each pore. This array must be
            *``Np``*-long, except when ``pores`` is also passed, where in that
            case they must be of equal length, with one value between 0 and 1,
            for each pore in the network. If a scalar is received it is applied
            to all pores. If nothing is passed, ``Pvals=1.0`` is assumed.

        Tvals : array_like, float
            The volume fraction of ``phase`` in each throat. This array must be
            *``Nt``*-long, except when ``throats`` is also passed, where in that
            case they must be of equal length, with one value between 0 and 1,
            for each throat in the network. If a scalar is received it is applied
            to all throats. If nothing is passed, ``Tvals=1.0`` is assumed.

        pores : array_like, int
            The location of pores for which the phase occupancy is to be set.

        throats : array_like, int
            The location of throats for which the phase occupancy is to be set.

        """
        Pvals = np.array(Pvals, ndmin=1)
        Tvals = np.array(Tvals, ndmin=1)
        pores =np.array(pores, ndmin=1)
        throats =np.array(throats, ndmin=1)

        # Check for size consistency of the arguments
        if Pvals.size and pores.size:
            if Pvals.size != pores.size:
                raise Exception("Pvals and pores must be the same size.")
        if Tvals.size and throats.size:
            if Tvals.size != throats.size:
                raise Exception("Tvals and throats must be the same size.")

        # Check if the passed phase is already part of MultiPhase object
        if phase not in self.project:
            raise Exception(f"{phase.name} doesn't belong to this project")
        # Add the passed phase to MultiPhase object if not found
        else:
            if phase.name not in self.settings['phases']:
                self.settings['phases'].append(phase.name)

        # Check for value consistency of the arguments
        if np.any(Pvals > 1.0) or np.any(Pvals < 0.0):
            logger.warning('Received Pvals contain volume fractions outside '
                           + 'the range of 0 to 1')
        if np.any(Tvals > 1.0) or np.any(Tvals < 0.0):
            logger.warning('Received Tvals contain volume fractions outside '
                           + 'the range of 0 to 1')

        if Pvals.size and not pores.size:
            pores = self.pores()
        if Tvals.size and not throats.size:
            throats = self.throats()

        if pores.size:
            Pvals = Pvals if Pvals.size else 1.0
            self['pore.occupancy.' + phase.name][pores] = Pvals
        if throats.size:
            Tvals = Tvals if Tvals.size else 1.0
            self['throat.occupancy.' + phase.name][throats] = Tvals

        self.regenerate_models(propnames=f"throat.occupancy.{phase.name}")
        self._update_occupancy()
