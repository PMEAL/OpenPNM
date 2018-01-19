from openpnm.core import Base, Workspace, logging, ModelsMixin
from openpnm.utils import PrintableDict
logger = logging.getLogger(__name__)
ws = Workspace()


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

    _prefix = 'phase'

    def __init__(self, network, name=None):
        super().__init__(name=name, simulation=network.simulation)
        self.settings = PrintableDict()
        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = network['pore.all']
        self['throat.all'] = network['throat.all']
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
            return self._interleave_data(key, self.simulation.physics.values())
        else:
            return super().__getitem__(key)
