import time
from openpnm.core import Workspace
ws = Workspace()


class Simulation(dict):

    def __init__(self, name=None):
        super().__init__()
        if name is None:
            name = 'sim'+str(len(ws.keys())+1).zfill(3)
        self.name = name
        self._grid = {}
        ws.update({name: self})

    def find_phase(self, physics):
        mro = [c.__name__ for c in physics.__class__.__mro__]
        if 'GenericPhase'in mro:
            return physics
        for g in self.geometries.values():
            for p in self.phases.values():
                if physics.name == self.grid[g.name][p.name]:
                    return p

    def find_geometry(self, physics):
        for g in self.geometries.values():
            for p in self.phases.values():
                if physics.name == self.grid[g.name][p.name]:
                    return g

    def find_physics(self, geometry, phase):
        name = self.grid[geometry.name][phase.name]
        return self[name]

    def add_network(self, network):
        if self.network is not None:
            raise Exception("This simulation already has a Network")
        self.update({network.name: network})

    def add_phase(self, phase):
        if phase in self.phases.values():
            raise Exception("A Phase with that name has already been added")
        self.update({phase.name: phase})

    def add_geometry(self, geometry):
        if geometry in self.geometries.values():
            raise Exception("A Geometry with that name has already been added")
        self.update({geometry.name: geometry})

    def add_physics(self, physics, geometry, phase):
        if physics in self.physics.values():
            raise Exception("A Physics with that name has already been added")
        self.update({physics.name: physics})

    def add_item(self, item):
        if item in self:
            raise Exception("An object with that name is already present " +
                            "in this simulation")

    def add_algorithm(self, algorithm):
        self.update({algorithm.name: algorithm})

    def validate_name(self, name):
        flag = True
        if name in self.keys():
            flag = False
        return flag

    def _get_net(self):
        for item in self.values():
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericNetwork' in mro:
                return item
    network = property(fget=_get_net)

    def _get_geoms(self):
        _dict = {}
        for item in self.values():
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericGeometry' in mro:
                _dict.update({item.name: item})
        return _dict
    geometries = property(fget=_get_geoms)

    def _get_phases(self):
        _dict = {}
        for item in self.values():
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericPhase' in mro:
                _dict.update({item.name: item})
        return _dict
    phases = property(fget=_get_phases)

    def _get_physics(self):
        _dict = {}
        for item in self.values():
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericPhysics' in mro:
                _dict.update({item.name: item})
        return _dict
    physics = property(fget=_get_physics)

    def _get_algorithms(self):
        _dict = {}
        for item in self.values():
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericAlgorithm' in mro:
                _dict.update({item.name: item})
        return _dict
    algorithms = property(fget=_get_algorithms)

    def _set_comments(self, string):
        if hasattr(self, '_comments') is False:
            self._comments = {}
        self._comments[time.strftime('%c')] = string

    def _get_comments(self):
        if hasattr(self, '_comments') is False:
            self._comments = {}
        for key in list(self._comments.keys()):
            print(key, ': ', self._comments[key])

    comments = property(fget=_get_comments, fset=_set_comments)

    def _get_grid(self):
        net = self.network
        grid = Grid()
        for geo in self.geometries.values():
            grid[geo.name] = {}
            Pg = net.pores(geo.name)[0]
            for phase in self.phases.values():
                grid[geo.name][phase.name] = ''
                for phys in self.physics.values():
                    if 'pore.'+phys.name in phase.keys():
                        if phase.pores(phys.name)[0] == Pg:
                            grid[geo.name][phase.name] = phys.name
        self._grid = grid
        return grid

    grid = property(fget=_get_grid)


class Grid(dict):

    def _get_sim(self):
        for sim in ws.values():
            if sim._grid is self:
                return sim

    def _get_geometries(self):
        sim = self._get_sim()
        return list(sim.geometries.keys())

    geometries = property(fget=_get_geometries)

    def _get_phases(self):
        sim = self._get_sim()
        return list(sim.phases.keys())

    phases = property(fget=_get_phases)

    def _get_net(self):
        sim = self._get_sim()
        return sim.network

    network = property(fget=_get_net)

    def __repr__(self):
        s = []
        hr = '―'*(16*(len(self.phases)+1))
        s.append(hr)
        fmt = ["| {"+str(i)+":^13} " for i in range(len(self.phases))]
        phases = [item for item in self.phases]
        s.append('| {0:^13}'.format(self.network.name) +
                 ''.join(fmt).format(*phases) + '|')
        s.append(hr)
        for geo in self.geometries:
            ind = '| {0:^13}'.format(geo)
            row = list(self[geo].values())
            s.append(ind + ''.join(fmt).format(*row) + '|')
            s.append(hr)
        return '\n'.join(s)
