from openpnm.core import Base, Workspace, ModelsMixin, logging
import scipy as sp
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericPhysics(Base, ModelsMixin):
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
        should be applied.

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that include the class name and a random string.

    """

    _prefix = 'phys'

    def __init__(self, network=None, phase=None, geometry=None, name=None):
        super().__init__(name=name, simulation=network.simulation)
        # Create a settings attribute
        self.settings['local_data'] = network.simulation.settings['local_data']
        # Initialize a label dictionary in the associated phase and network
        phase['pore.'+self.name] = False
        phase['pore.'+self.name][network.pores(geometry.name)] = True
        phase['throat.'+self.name] = False
        phase['throat.'+self.name][network.throats(geometry.name)] = True
        self.update({'pore.all':
                     sp.ones(shape=sp.size(geometry.Ps), dtype=bool)})
        self.update({'throat.all':
                     sp.ones(shape=sp.size(geometry.Ts), dtype=bool)})

    def __getitem__(self, key):
        element = key.split('.')[0]
        if key.split('.')[-1] == '_id':
            phase = self.simulation.find_phase(self)
            inds = phase._get_indices(element=element, labels=self.name)
            vals = phase[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] in [self.name]:
            vals = self[element+'.all']
        # Get prop or label if present
        elif key in self.keys():
            vals = super(Base, self).__getitem__(key)
        # Otherwise retrieve from network
        else:
            phase = self.simulation.find_phase(self)
            inds = phase._get_indices(element=element, labels=self.name)
            vals = phase[key][inds]
        return vals

    def __setitem__(self, key, value):
        if self.settings['local_data']:
            super().__setitem__(key, value)
        else:
            phase = self.simulation.find_phase(self)
            element = self._parse_element(key.split('.')[0], single=True)
            inds = self._map(ids=self[element+'._id'], element=element,
                             filtered=True)
            # If array not in network, create it first
            if key not in phase.keys():
                if value.dtype == bool:
                    phase[key] = False
                else:
                    phase[key] = sp.zeros_like(value)*sp.nan
            phase[key][inds] = value
