from openpnm.phases import GenericPhase as GenericPhase
import openpnm.models.phases as _models
import numpy as _np
from openpnm.utils import logging
logger = logging.getLogger(__name__)

def_sets = {'pore_occupancy': 'pore.occupancy',
            'throat_occupancy': 'throat.occupancy',
            'phases': []
            }


class Multiphase(GenericPhase):
    r"""
    Creates Phase object that represents a multiphase system consisting of
    a given list of OpenPNM Phase objects.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    project : OpenPNM Project object, optional
        The Project with which this phase should be associted.  If a
        ``network`` is given then this is ignored and the Network's project is
        used.  If a ``network`` is not given then this is mandatory.

    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.  If
        no name is given, one is generated.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> air = op.phases.Air(network=pn)
    >>> water = op.phases.Water(network=pn)
    >>> multiphase = op.phases.Multiphase(network=pn, phases=[air, water])

    """
    def __init__(self, phases=[], settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(def_sets)
        self.settings.update(settings)

        # Initialize occupancy 'dummy' arrays to show up in props
        self.update({self.settings['pore_occupancy']:
                     _np.nan*_np.ones([self.Np, ])})
        self.update({self.settings['throat_occupancy']:
                     _np.nan*_np.ones([self.Nt, ])})

        # Add any supplied phases to the phases list
        for phase in phases:
            self.add_phase(phase)

    def __getitem__(self, key):
        # If occupancy is requested, build it from component phases
        if key in [self.settings['pore_occupancy'],
                   self.settings['throat_occupancy']]:
            element = key.split('.')[0]
            vals = _np.zeros_like(self._get_indices(element=element),
                                  dtype=float)
            for phase in self.phases:
                vals = vals + phase[key].astype(float)
            return vals
        # Otherwise just retrieve requested array
        vals = super().__getitem__(key)
        return vals

    def __setitem__(self, key, value):
        # Prevent attempts to write to occupancy arrays
        if key in [self.settings['pore_occupancy'],
                   self.settings['throat_occupancy']]:
            raise Exception("Cannot set occupancy on MultiPhase, it must be " +
                            "set on component phases")
        # Let all others through
        super().__setitem__(key, value)

    def add_phase(self, phase):
        r"""
        Adds a given phase to the list of phases composing this mixture

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to add to the list
        """
        if phase not in self.project:
            raise Exception(f"{phase.name} doesn't belong to this project")
        if phase in self.phases:
            logger.warning(phase.name + ' already present')
        else:
            self.settings['phases'].append(phase.name)
        if self.settings['pore_occupancy'] not in phase.keys():
            phase[self.settings['pore_occupancy']] = 0.0
        if self.settings['throat_occupancy'] not in phase.keys():
            phase[self.settings['throat_occupancy']] = 0.0

    def _get_phases(self):
        phases = [self.project[item] for item in self.settings['phases']]
        return phases

    phases = property(fget=_get_phases)

    def interleave_data(self, prop, mixing_rule='fractional'):
        r"""
        Gathers property values from component phases to build a single array

        Parameters
        ----------
        prop : string
            The property to be retrieved.
        mixing_rule : string
            Indicates how properties from different phase should be blended.
            Options are:

            'fractional' - Multiplies the ``prop`` of each phase by the
            fractional occupancy of that phase to produce a weighted average.
            This also works with boolean values to apply a mask.

            'custom' - Multiplies the ``prop`` by the specified ``propname``
            to produce a weighted average, such as molar mass or density.
            Can be any custom array.

        Returns
        -------
        array : ND-array
            An array containing the specified property retrieved from each
            component phase and assembled based on the specified mixing rule

        """
        element = prop.split('.')[0]
        vals = _np.zeros([self._count(element=element)], dtype=float)
        for phase in self.phases:
            vals += phase[prop]*phase[self.settings[element + '_occupancy']]
        return vals

    # Models of constituent phases must update before those of Multiphase
    def regenerate_models(self, **kwargs):
        r"""
        Regenerate models associated with the Multiphase object

        This method works by first regenerating the models associated with the
        constituent phases, and then regenerating Multiphase models.

        """
        # Regenerate all phases within mixture
        for phase in self.phases:
            phase.regenerate_models(**kwargs)
        # Regenerate mixture
        super().regenerate_models(self, **kwargs)

    def set_occupancy(self, phase, pores=[], throats=[], values=1.0):
        r"""
        Specify occupancy of a phase in each pore and/or throat

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase whose occupancy is being specified

        pores : array_like
            The pore indices where occupancy is being set

        throats : array_like
            The thraot indices where occupancy is being set

        values : array_like
            The occupancy to apply to the given pores and/or throats.  Note
            that if a scalar value is given it will be applied to *both* the
            specified pores and throats (if given).  If an array is given it
            must match the length of the given pore or throat indices.

        Notes
        -----
        This method actually only writes the given pore/throat occupancy on
        the ``phase`` object itself, not the Multiphase object.  It is
        identical to setting ``phase['pore.occupancy'][pores] = values``.

        """
        values = _np.array(values)
        pores = _np.array(pores)
        throats = _np.array(throats)
        # Ensure phase is part of multiphase
        if phase.name not in self.settings['phases']:
            raise Exception('Provided Phase is not part of this MultiPhase')
        if pores.size > 0:
            phase[self.settings['pore_occupancy']][pores] = values
        if throats.size > 0:
            phase[self.settings['pore_occupancy']][throats] = values
