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

    def __init__(self, network=None, project=None, pores=[], throats=[],
                 settings={}, **kwargs):
        # Define some default settings
        self.settings.update({'prefix': 'geo'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)

        # Deal with network or project arguments
        if network is not None:
            project = network.project

        super().__init__(project=project, **kwargs)

        self.add_locations(pores=pores, throats=throats)

    def __getitem__(self, key):
        net = self.project.network
        element = key.split('.')[0]
        # Get uuid from network
        if key.split('.')[-1] == '_id':
            inds = net._get_indices(element=element, labels=self.name)
            vals = net[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] in [self.name]:
            vals = self[element+'.all']
        # Apply logic in the __missing__ method
        elif key not in self.keys():
            vals = self.__missing__(key)
        else:
            vals = super(Base, self).__getitem__(key)
        return vals

    def __missing__(self, key):
        net = self.project.network
        element = key.split('.')[0]
        # If key not available try running model
        if key in self.models.keys():
            print('GenericGeometry: ' + key + ' missing, running model')
            self.regenerate_models(propnames=[key])
            vals = self.__getitem__(key)
        # If not found on network a key error will be raised
        else:
            inds = net._get_indices(element=element, labels=self.name)
            vals = super(Base, net).__getitem__(key)[key][inds]
        return vals

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
        network = self.project.network
        indices = network._parse_indices(indices)

        net = self.project.network
        proj = self.project

        if element+'.'+self.name not in net.keys():
            net[element+'.'+self.name] = False

        if mode == 'add':
            # Ensure indices are not already assigned to another object
            temp = sp.zeros(shape=[net._count(element=element), ], dtype=bool)
            for item in proj.geometries().keys():
                temp += net[element+'.'+item]
            if sp.any(temp[indices]):
                raise Exception('Some of the given '+element+' are already ' +
                                'assigned to an existing object')
            set_flag = True
        elif mode == 'drop':
            set_flag = False

        # Change lables of all associated physics in their respective phases
        for phase in proj.phases().values():
            phys = proj.find_physics(geometry=self, phase=phase)
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
        return self.project.network
