import scipy as sp
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
    >>> Ps = pn.pores()  # Get all pores
    >>> Ts = pn.throats()  # Get all throats
    >>> geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    """

    def __init__(self, network, pores=[], throats=[], name=None):
        super().__init__(name=name, simulation=network.simulation)
        logger.name = self.name

        # Initialize a label dictionary in the associated network
        network['pore.'+self.name] = False
        network['throat.'+self.name] = False
        self['pore.all'] = sp.ones(shape=sp.size(pores), dtype=bool)
        self['throat.all'] = sp.ones(shape=sp.size(throats), dtype=bool)
        network['pore.'+self.name] = False
        network['pore.'+self.name][pores] = True
        network['throat.'+self.name] = False
        network['throat.'+self.name][throats] = True

    def __getitem__(self, key):
        net = self.simulation.network
        element = key.split('.')[0]
        # Convert self.name into 'all'
        if key.split('.')[-1] == self.name:
            key = element + '.all'
        if key in list(self.keys()):  # Look for data on self...
            return super(GenericGeometry, self).__getitem__(key)
        if key == 'throat.conns':  # Handle specifically
            [P1, P2] = net['throat.conns'][net[element+'.'+self.name]].T
            Pmap = sp.zeros((self._net.Np,), dtype=int) - 1
            Pmap[self._net.pores(self.name)] = self.Ps
            conns = sp.array([Pmap[P1], Pmap[P2]]).T
            # Replace -1's with nans
            if sp.any(conns == -1):
                conns = sp.array(conns, dtype=object)
                conns[sp.where(conns == -1)] = sp.nan
            return conns
        else:  # ...Then check Network
            return net[key][net[element+'.'+self.name]]

    def set_locations(self, pores=None, throats=None, mode='add'):
        r"""
        """
        net = self.simulation.network
        sim = self.simulation
        if pores is not None:
            element = 'pore'
            indices = self._parse_indices(pores)
        elif throats is not None:
            element = 'throats'
            indices = self._parse_indices(throats)
        else:
            raise Exception('Can\'t set pores and throats at the same time')

        if mode == 'add':
            # Ensure indices are not already assigned to another object
            temp = sp.zeros(shape=[net._count(element=element), ], dtype=bool)
            for item in sim.geometries.keys():
                temp += net[element+'.'+item]
            if sp.any(temp[indices]):
                raise Exception('Some of the given '+element+' are already ' +
                                'assigned to an existing object')

            # Create new 'all' label for new size
            new_len = self._count(element=element) + sp.size(indices)
            self.update({element+'.all': sp.ones((new_len, ), dtype=bool)})

            # Update indices in network and ph
            inds_orig = net._get_indices(element=element, labels=self.name)
            if element+'.'+self.name not in net.keys():
                net[element+'.'+self.name] = False
            net[element+'.'+self.name][indices] = True
            inds_new = net._get_indices(element=element, labels=self.name)

            # Increase size of labels (add False at new indices)
            labels = self.labels()
            labels.remove(element+'.all')
            for item in labels:
                if item.split('.')[0] == element:
                    net[element+'.'+'blank'] = False
                    net[element+'.'+'blank'][inds_orig] = self[item]
                    self[item] = net[element+'.'+'blank'][inds_new]
            net.pop(element+'.'+'blank', None)

        if mode == 'remove':
            # Change the labeling in the boss object
            net[element+'.'+self.name][indices] = False
            # Convert network indices to obj-specific indices
            obj_inds = net._map(element=element,
                                locations=indices,
                                target=self)
            keep = ~self._tomask(indices=obj_inds, element=element)
            for item in list(self.keys()):
                if item.split('.')[0] == element:
                    temp = self[item][keep]
                    self.update({item: temp})
