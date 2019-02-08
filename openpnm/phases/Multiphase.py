from openpnm.phases import GenericPhase as _GenericPhase
import openpnm.models.phases as _models
import numpy as _np


class Multiphase(_GenericPhase):
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

    References
    ----------
    The pore scale models for this class are taken from its constituent phases.

    """
    def __init__(self, phases, occupancy="pore.occupancy", **kwargs):
        super().__init__(**kwargs)
        self.phases = phases
        self.occupancy = occupancy

        # Ensure all phases belong to a single network
        for phase in self.phases:
            if phase.project.network is not self.project.network:
                raise Exception(f"{phase.name} doesn't belong to this network")

        # Enforce zero occupancy for all phases if not present
        for phase in self.phases:
            if occupancy not in phase.props():
                phase[occupancy] = False

        # Find shared props
        props = self.phases[0].props()
        for phase in self.phases:
            props = list(set(props).intersection(phase.props()))

        # By default, the mixing rule is based on volume occupancy
        for prop in props:
            self.add_model(propname=prop, model=_models.misc.mix_and_match,
                           prop=prop+"_", occupancy=occupancy, phases=self.phases)

    # Models of constituent phases must update before those of Multiphase
    def regenerate_models(self, **kwargs):
        r"""
        Regenerate models associated with the Multiphase object.

        This method works by first regenerating the models associated with the
        consituent phases, and then regenerating Multiphase models.

        """
        # Regenerate all phases within Mixture
        for phase in self.phases:
            phase.regenerate_models(*kwargs)
        # Regenerate Mixture phase
        super().regenerate_models(self, *kwargs)

    def set_occupancy(self, phase, pores, values):
        r"""

        """
        values = _np.array(values, dtype=bool)
        phase[self.occupancy][pores] = values
        other_phases = [elem for elem in self.phases if elem is not phase]
        for elem in other_phases:
            elem[self.occupancy][pores] = ~values
