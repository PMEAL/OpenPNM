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
    >>> print(m)
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    openpnm.phases.Multiphase : multi
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    #     Properties                                    Valid Values
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1     pore.occupancy.all                              125 / 125
    2     throat.occupancy.all                            300 / 300
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    #     Labels                                        Assigned Locations
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    1     pore.all                                      125
    2     pore.occupancy.air                            75
    3     pore.occupancy.water                          50
    4     throat.all                                    300
    5     throat.occupancy.air                          195
    6     throat.occupancy.water                        105
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    Air and water have uniform viscosity values throughout, so the peak-to-peak
    distances is 0, while the mixture phase has the viscosity of air in some
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
        self.settings.update({'phases': [],
                              })
        self.settings.update(settings)

        self['pore.occupancy.all'] = np.zeros(self.Np, dtype=float)
        self['throat.occupancy.all'] = np.zeros(self.Nt, dtype=float)

        self.pop('pore.temperature', None)
        self.pop('pore.pressure', None)

        # Add any supplied phases to the phases list
        for phase in phases:
            self.settings['phases'].append(phase.name)

    def __getitem__(self, key):
        try:
            vals = super().__getitem__(key)
        except KeyError:
            vals = self.interleave_data(key)
        return vals

    def _update_occupancy(self):
        # Update occupancy.all
        self['pore.occupancy.all'] = 0.0
        dict_ = list(self['pore.occupancy'].values())
        if len(dict_) > 1:
            self['pore.occupancy.all'] = np.sum(dict_, axis=0)
        self['throat.occupancy.all'] = 0.0
        dict_ = list(self['throat.occupancy'].values())
        if len(dict_) > 1:
            self['throat.occupancy.all'] = np.sum(dict_, axis=0)

    def _get_phases(self):
        phases = [self.project[item] for item in self.settings['phases']]
        return phases

    phases = property(fget=_get_phases)

    def interleave_data(self, prop):
        r"""
        Gathers property values from component phases to build a single array

        If the requested ``prop`` is not on this MultiPhase, then a search is
        conducted on all associated physics objects, and values from each
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
        element = prop.split('.')[0]
        if np.any(self[element + '.occupancy.all'] != 1.0):
            self._update_occupancy()
            if np.any(self[element + '.occupancy.all'] != 1.0):
                raise Exception('Occupancy does not add to unity in all ' +
                                element + 's')
        vals = np.zeros([self._count(element=element)], dtype=float)
        try:
            for phase in self.phases:
                vals += phase[prop]*self[element + '.occupancy.' + phase.name]
        except KeyError:
            vals = super().interleave_data(prop)
        return vals

#    # Models of constituent phases must update before those of Multiphase
#    def regenerate_models(self, **kwargs):
#        r"""
#        Regenerate models associated with the Multiphase object
#
#        This method works by first regenerating the models associated with the
#        constituent phases, and then regenerating Multiphase models.
#
#        """
#        # Regenerate all phases within mixture
#        for phase in self.phases:
#            phase.regenerate_models(**kwargs)
#        # Regenerate mixture
#        super().regenerate_models(self, **kwargs)

    def set_occupancy(self, phase, Pvals=[], Tvals=[]):
        r"""
        Specify occupancy of a phase in each pore and/or throat

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase whose occupancy is being specified

        Pvals : array_like
            The volume fraction of ``phase`` in each pore.  This array must
            be *Np*-long, with one value between 0 and 1, for each pore in the
            network.  If a scalar is received it is applied to all pores.

        Tvals : array_like
            The volume fraction of ``phase`` in each throat.  This array must
            be *Nt*-long, with one value between 0 and 1, for each throat in
            the network.  If a scalar is received it is applied to all throats.

        """
        Pvals = np.array(Pvals, ndmin=1)
        Tvals = np.array(Tvals, ndmin=1)
        if phase not in self.project:
            raise Exception(f"{phase.name} doesn't belong to this project")
        else:
            if phase.name not in self.settings['phases']:
                self.settings['phases'].append(phase.name)
        if np.any(Pvals > 1.0) or np.any(Pvals < 0.0):
            logger.warning('Received Pvals contain volume fractions outside ' +
                           'the range of 0 to 1')
        if np.any(Tvals > 1.0) or np.any(Tvals < 0.0):
            logger.warning('Received Tvals contain volume fractions outside ' +
                           'the range of 0 to 1')
        if Pvals.size:
            self['pore.occupancy.' + phase.name] = Pvals
        if Tvals.size:
            self['throat.occupancy.' + phase.name] = Tvals
        self._update_occupancy()
