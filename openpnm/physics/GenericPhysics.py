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

    def __init__(self, network=None, phase=None, geometry=None, name=None):
        super().__init__(name=name, simulation=network.simulation)

        # Initialize a label dictionary in the associated phase and network
        self['pore.all'] = sp.ones(shape=sp.size(geometry.Ps), dtype=bool)
        self['throat.all'] = sp.ones(shape=sp.size(geometry.Ts), dtype=bool)
        phase['pore.'+self.name] = False
        phase['pore.'+self.name][network.pores(geometry.name)] = True
        phase['throat.'+self.name] = False
        phase['throat.'+self.name][network.throats(geometry.name)] = True

    def __getitem__(self, key):
        phase = self.simulation.find_phase(self)
        element = key.split('.')[0]
        if key.split('.')[-1] == '_id':
            inds = phase._get_indices(element=element, labels=self.name)
            return phase[element+'._id'][inds]
        # Convert self.name into 'all'
        if key.split('.')[-1] == self.name:
            key = element + '.all'
        if key in self.keys():  # Look for data on self...
            return super(GenericPhysics, self).__getitem__(key)
        else:  # ...Then check Phase
            return phase[key][phase[element + '.' + self.name]]
