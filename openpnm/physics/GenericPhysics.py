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

    def __init__(self, network, phase, geometry, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'phys')
        self.settings.update({'boss': phase.name})
        self.settings.update(settings)
        super().__init__(Np=geometry.Np, Nt=geometry.Nt,
                         simulation=network.simulation, **kwargs)
        self.settings['local_data'] = self.simulation.settings['local_data']
        # Initialize a label array in the associated phase
        phase['pore.'+self.name] = False
        phase['pore.'+self.name][network.pores(geometry.name)] = True
        phase['throat.'+self.name] = False
        phase['throat.'+self.name][network.throats(geometry.name)] = True

    def __getitem__(self, key):
        element = key.split('.')[0]
        boss = self.simulation[self.settings['boss']]
        if key.split('.')[-1] == '_id':
            inds = boss._get_indices(element=element, labels=self.name)
            vals = boss[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] in [self.name]:
            vals = self[element+'.all']
        # Get prop or label if present
        elif key in self.keys():
            vals = super(Base, self).__getitem__(key)
        # Otherwise retrieve from network
        else:
            inds = boss._get_indices(element=element, labels=self.name)
            vals = boss[key][inds]
        return vals

    def __setitem__(self, key, value):
        if self.settings['local_data']:
            super().__setitem__(key, value)
        else:
            boss = self.simulation[self.settings['boss']]
            element = self._parse_element(key.split('.')[0], single=True)
            inds = boss._map(ids=self[element+'._id'], element=element,
                             filtered=True)
            # If array not in phase, create it first
            if key not in boss.keys():
                if value.dtype == bool:
                    boss[key] = False
                else:
                    dtype = value.dtype
                    if dtype.name == 'object':
                        boss[key] = sp.zeros(1, dtype=object)
                    else:
                        Nt = len(boss[element+'.all'])
                        dim = sp.size(value[0])
                        if dim > 1:
                            arr = sp.zeros(dim, dtype=dtype)
                            temp = sp.tile(arr, reps=(Nt, 1))*sp.nan
                        else:
                            temp = sp.zeros(Nt)*sp.nan
                        boss[key] = temp
            boss[key][inds] = value
