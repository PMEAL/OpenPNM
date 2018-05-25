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

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        # Define some default settings
        self.settings.update({'prefix': 'phase'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project

        super().__init__(project=project, **kwargs)

        # If project has a network object, adjust pore and throat sizes
        network = self.project.network
        if network:
            self['pore.all'] = ones((network.Np, ), dtype=bool)
            self['throat.all'] = ones((network.Nt, ), dtype=bool)

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __setitem__(self, key, value):
        if self.project:
            for item in self.project.find_physics(phase=self):
                exclude = {'pore.all', 'throat.all'}
                if key in set(item.keys()).difference(exclude):
                    raise Exception(key+' already exists on '+item.name)
        super().__setitem__(key, value)

    def __getitem__(self, key):
        element = key.split('.')[0]
        # Deal with special keys first
        if key.split('.')[-1] == '_id':
            net = self.project.network
            return net[element+'._id']
        if key.split('.')[-1] == self.name:
            return self[element+'.all']
        # Now get values if present, or regenerate them
        vals = self.get(key)
        if vals is None:
            physics = self.project.find_physics(phase=self)
            vals = self._interleave_data(key, physics)
        return vals
