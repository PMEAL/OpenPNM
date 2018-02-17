from openpnm.core import Base, Workspace, logging, ModelsMixin
from openpnm.utils import PrintableDict
logger = logging.getLogger(__name__)
ws = Workspace()
from numpy import ones


class GenericPhase(Base, ModelsMixin):
    r"""
    Base class to generate a generic phase object.  The user must specify
    models and parameters for all the properties they require. Classes for
    several common phases are included with OpenPNM and can be found under
    ``openpnm.phases``.

    Parameters
    ----------
    network : openpnm Network object
        The network to which this Phase should be attached

    name : str, optional
        A unique string name to identify the Phase object, typically same as
        instance name but can be anything.

    """

    def __init__(self, network=None, simulation=None, settings={}, **kwargs):
        # Define some default settings
        self.settings.update({'prefix': 'phase'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or simulation arguments
        if network is not None:
            simulation = network.simulation

        super().__init__(simulation=simulation, **kwargs)

        # If simulation has a network object, adjust pore and throat sizes
        if simulation.network:
            self['pore.all'] = ones((simulation.network.Np, ), dtype=bool)
            self['throat.all'] = ones((simulation.network.Nt, ), dtype=bool)

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __getitem__(self, key):
        element = key.split('.')[0]
        if key.split('.')[-1] == '_id':
            net = self.simulation.network
            return net[element+'._id']
        if key.split('.')[-1] == self.name:
            return self[element+'.all']
        if key not in self.keys():
            logger.debug(key + ' not on Phase, constructing data from Physics')
            names = self.simulation.find_physics(phase=self)
            physics = [self.simulation.physics()[i] for i in names]
            return self._interleave_data(key, physics)
        else:
            return super().__getitem__(key)
