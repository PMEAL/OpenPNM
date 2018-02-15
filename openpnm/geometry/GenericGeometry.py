import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from openpnm.core import Base, ModelsMixin
from openpnm.core import logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericGeometry(Base, ModelsMixin):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object to which this Geometry applies.

    pores and/or throats : array_like
        The list of pores and throats where this Geometry applies.

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this Geometry applies.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> Ps = pn.pores('all')  # Get all pores
    >>> Ts = pn.throats('all')  # Get all throats
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    """

    def __init__(self, network, pores=[], throats=[], settings={}, **kwargs):
        self.settings.setdefault('prefix', 'geo')
        self.settings.update(settings)
        super().__init__(simulation=network.simulation, **kwargs)
        self.settings['local_data'] = self.simulation.settings['local_data']
        self.add_locations(pores=pores, throats=throats)

    def __getitem__(self, key):
        net = self.simulation.network
        element = key.split('.')[0]
        inds = net._get_indices(element=element, labels=self.name)
        # Get uuid from network
        if key.split('.')[-1] == '_id':
            vals = net[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] in [self.name]:
            vals = self[element+'.all']
        # Get prop or label if present
        elif key in self.keys():
            vals = super(Base, self).__getitem__(key)
        else:
            # If not found on network a key error will be raised
            vals = net[key][inds]
        return vals

    def __setitem__(self, key, value):
        if self.settings['local_data']:
            super().__setitem__(key, value)
        else:
            network = self.simulation.network
            element = self._parse_element(key.split('.')[0], single=True)
            inds = network._map(ids=self[element+'._id'], element=element,
                                filtered=True)
            # If array not in network, create it first
            if key not in network.keys():
                if value.dtype == bool:
                    network[key] = False
                else:
                    dtype = value.dtype
                    if dtype.name == 'object':
                        network[key] = np.zeros(1, dtype=object)
                    else:
                        Nt = len(network[element+'.all'])
                        dim = np.size(value[0])
                        if dim > 1:
                            arr = np.zeros(dim, dtype=dtype)
                            temp = np.tile(arr, reps=(Nt, 1))*np.nan
                        else:
                            temp = np.zeros(Nt)*np.nan
                        network[key] = temp

            network[key][inds] = value

    def add_locations(self, pores=[], throats=[]):
        r"""
        """
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='add')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='add')

    def drop_locations(self, pores=[], throats=[]):
        r"""
        """
        if len(pores) > 0:
            self._set_locations(element='pore', indices=pores, mode='drop')
        if len(throats) > 0:
            self._set_locations(element='throat', indices=throats, mode='drop')

    def _set_locations(self, element, indices=[], mode='add'):
        r"""
        """
        element = self._parse_element(element=element, single=True)
        # Use the network's _parse_indices, since indicies could be 'network'
        # length boolean masks
        network = self.simulation.network
        indices = network._parse_indices(indices)

        net = self.simulation.network
        sim = self.simulation

        if element+'.'+self.name not in net.keys():
            net[element+'.'+self.name] = False

        if mode == 'add':
            # Ensure indices are not already assigned to another object
            temp = sp.zeros(shape=[net._count(element=element), ], dtype=bool)
            for item in sim.geometries().keys():
                temp += net[element+'.'+item]
            if sp.any(temp[indices]):
                raise Exception('Some of the given '+element+' are already ' +
                                'assigned to an existing object')
            set_flag = True
        elif mode == 'drop':
            set_flag = False

        # Change lables of all associated physics in their respective phases
        for phase in sim.phases().values():
            phys = sim.find_physics(geometry=self, phase=phase)
            if phys:
                if element+'.'+phys.name not in phase.keys():
                    phase[element+'.'+phys.name] = False
                phase[element+'.'+phys.name][indices] = set_flag
                temp = sp.ones(shape=(sp.sum(phase[element+'.'+phys.name]),),
                               dtype=bool)
                phys.update({element+'.all': temp})

        # Change labels in the network
        net[element+'.'+self.name][indices] = set_flag
        temp = sp.ones(shape=(sp.sum(net[element+'.'+self.name]),), dtype=bool)
        self.update({element+'.all': temp})
        inds = net._get_indices(element=element, labels=self.name)
        for item in self.props(element=element):
            self.update({item: net[item][inds]})

    def show_hist(self, props=['pore.diameter'], bins=20, fig=None, **kwargs):
        if fig is None:
            fig = plt.figure()
        if type(props) is str:
            props = [props]
        N = len(props)
        if N == 1:
            r = 1
            c = 1
        elif N < 4:
            r = 1
            c = N
        else:
            r = int(sp.ceil(N**0.5))
            c = int(sp.floor(N**0.5))

        for i in range(len(props)):
            plt.subplot(r, c, i+1)
            plt.hist(self[props[i]], bins=bins, **kwargs)

    @property
    def network(self):
        r"""
        A shortcut to get a handle to the associated network
        There can only be one so this works
        """
        return self.simulation.network
