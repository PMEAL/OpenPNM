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

        if network is not None:
            network['pore.'+self.name] = False
            network['throat.'+self.name] = False

        self.add_locations(pores=pores, throats=throats)

    def __getitem__(self, key):
        # Find boss object (either phase or network)
        element = key.split('.')[0]
        if self._isa('phase'):
            boss = self.project.find_phase(self)
        else:
            boss = self.project.network
        # Deal with a few special key items
        if key.split('.')[-1] == '_id':
            inds = boss._get_indices(element=element, labels=self.name)
            return boss[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] in [self.name]:
            return self[element+'.all']
        # Now get values if present, or regenerate them
        vals = self.get(key)
        if vals is None:
            inds = boss._get_indices(element=element, labels=self.name)
            vals = boss[key][inds]
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
        if self._isa('physics'):
            boss = self.project.find_phase(physics=self)
        if self._isa('geometry'):
            boss = self.project.network
        element = self._parse_element(element=element, single=True)
        # Use the network's _parse_indices, since indicies could be 'network'
        # length boolean masks
        indices = boss._parse_indices(indices)

        # Add self's label to network if not present, to prevent errors below
        if element+'.'+self.name not in boss.keys():
            boss[element+'.'+self.name] = False

        # Find mask of existing locations (network indexing)
        mask = boss[element+'.'+self.name]
        # Update mask with new locations (either add or remove)
        if mode == 'add':
            mask = mask + boss._tomask(indices=indices, element=element)
        elif mode == 'drop':
            mask = mask ^ (boss._tomask(indices=indices, element=element))
        # Change size of all arrays on self
        for item in self.keys(element=element, mode='all'):
            self.update({item: boss[item][mask]})
        # Update label array in network
        boss[element+'.'+self.name] = mask

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
