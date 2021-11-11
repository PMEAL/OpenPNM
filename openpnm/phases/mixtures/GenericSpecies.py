from openpnm.phases import GenericPhase as GenericPhase
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class GenericSpecies(GenericPhase):
    r"""
    Creates Phase object that represents a single species in a multicomponent
    mixture system.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.
    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    """
    @property
    def mixture(self):
        for item in self.project:
            if hasattr(item, 'components'):
                for comp in item.components.values():
                    if self is comp:
                        return item
