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
        self.update({'pore.all': sp.ones(shape=len(pores), dtype=bool)})
        self.update({'throat.all': sp.ones(shape=len(throats), dtype=bool)})
        network['pore.'+self.name] = False
        network['pore.'+self.name][pores] = True
        network['throat.'+self.name] = False
        network['throat.'+self.name][throats] = True

    def __getitem__(self, key):
        net = self.simulation.network
        element = key.split('.')[0]
        inds = net._get_indices(element=element, labels=self.name)
        # Get uuid from network
        if key.split('.')[-1] == '_id':
            vals = net[element+'._id'][inds]
        # Convert self.name into 'all'
        elif key.split('.')[-1] == self.name:
            vals = self[element+'.all']
        # Get prop or label is present
        elif key in self.keys():
            vals = self[key]
        # Otherwise retrieve from network
        else:
            vals = net[key][inds]
        return vals

    def add_locations(self, pores=None, throats=None):
        r"""
        """
        if pores is not None:
            element = 'pore'
            indices = self._parse_indices(pores)
        elif throats is not None:
            element = 'throat'
            indices = self._parse_indices(throats)
        else:
            raise Exception('Can\'t set pores and throats at the same time')

        net = self.simulation.network
        sim = self.simulation
        physics = []
        for phase in sim.phases.values():
            physics.append(sim.find_physics(geometry=self, phase=phase))

        orig_ids = self.map_pores(ids=sp.array(indices))

        # Ensure indices are not already assigned to another object
        temp = sp.zeros(shape=[net._count(element=element), ], dtype=bool)
        for item in sim.geometries.keys():
            temp += net[element+'.'+item]
        if sp.any(temp[indices]):
            raise Exception('Some of the given '+element+' are already ' +
                            'assigned to an existing object')

        # Update indices in network and phases
        if element+'.'+self.name not in net.keys():
            net[element+'.'+self.name] = False
        net[element+'.'+self.name][indices] = True
        for phys in filter(None, physics):
            phase = sim.find_phase(phys)
            if element+'.'+phys.name not in phase.keys():
                phase[element+'.'+phys.name] = False
            phase[element+'.'+phys.name][indices] = True

        # Increase the size of all on geometry and associated physics objects
        new_len = self._count(element=element) + sp.size(indices)
        new_all = sp.ones((new_len, ), dtype=bool)
        self.update({element+'.all': new_all})
        for phys in filter(None, physics):
            phys.update({element+'.all': new_all})

        # Increase the size of prop and label arrays on geom and all physics
        # First make list of constant props and labels that need changing
        keys = self.props(element=element, mode='constants')
        keys = keys + self.labels(element=element)
        keys.pop(keys.index(element+'.all'))
        inds = self.map_pores(ids=orig_ids)
        for item in keys:
            logger.warning(item + ' now contains invalid values')
            if self[item].dtype == bool:
                sp.insert(self[item], inds, False, axis=0)
            else:
                sp.insert(self[item].astype(float), inds, sp.nan, axis=0)

        self.clear(mode='models', element=element)

    def drop_locations(self, pores=None, throats=None):
        r"""
        """
        if pores is not None:
            element = 'pore'
            indices = self._parse_indices(pores)
        elif throats is not None:
            element = 'throat'
            indices = self._parse_indices(throats)
        else:
            raise Exception('Can\'t set pores and throats at the same time')

        net = self.simulation.network
        sim = self.simulation
        physics = []
        for phase in sim.phases.values():
            physics.append(sim.find_physics(geometry=self, phase=phase))

        # Ensure given indices are actually assigned to current geometry
        inds = self._map(ids=net[element+'._id'][indices], element=element,
                         filtered=True)
        if len(inds) < len(indices):
            raise Exception('Some provided locations are not assigned to ' +
                            'this geometry')

        inds = self._tomask(indices=inds, element=element)
        # Remove values from prop and label arrays on geom and all physics
        keys = self.props(element=element) + self.labels(element=element)
        for item in keys:
            temp = self.pop(item)
            self.update({item: temp[~inds]})
        for phys in filter(None, physics):
                keys = phys.props(element=element) + \
                       phys.labels(element=element)
                for item in keys:
                    temp = phys.pop(item)
                    phys.update({item: temp[~inds]})

        # Change the labeling in the network and phases
        net[element+'.'+self.name][indices] = False
        for phase in sim.phases.values():
            phys = sim.find_physics(geometry=self, phase=phase)
            if phys:
                phase[element+'.'+phys.name][indices] = False
