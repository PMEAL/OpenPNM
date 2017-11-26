import time
from openpnm.core import Workspace
ws = Workspace()


class Simulation(list):

    def __init__(self, name=None):
        super().__init__()
        self._name = None
        self.name = name
        self._grid = {}
        ws.update({self.name: self})

    def _set_name(self, name):
        if name is None:
            name = ws._gen_name()
        if name in ws.keys():
            raise Exception("A simulation with that name already exists")
        if self.name is not None:
            old_name = self.name
            ws.pop(old_name, None)
            ws.update({name: self})
        self._name = name

    def _get_name(self):
        return self._name

    name = property(fget=_get_name, fset=_set_name)

    def __getitem__(self, key):
        for obj in self:
            if obj.name == key:
                return obj

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

    def validate_name(self, name):
        flag = True
        names = [i.name for i in self]
        if name in names:
            flag = False
        return flag

    def _get_net(self):
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericNetwork' in mro:
                return item

    network = property(fget=_get_net)

    def _get_geoms(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericGeometry' in mro:
                _dict.update({item.name: item})
        return _dict

    geometries = property(fget=_get_geoms)

    def _get_phases(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericPhase' in mro:
                _dict.update({item.name: item})
        return _dict

    phases = property(fget=_get_phases)

    def _get_physics(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericPhysics' in mro:
                _dict.update({item.name: item})
        return _dict

    physics = property(fget=_get_physics)

    def _get_algorithms(self):
        _dict = {}
        for item in self:
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

    def __str__(self):
        s = []
        hr = '―'*80
        s.append(hr)
        s.append(' {0:<15} '.format('Object Name') +
                 '{0:<65}'.format('Object ID'))
        s.append(hr)
        for item in self:
            s.append(' {0:<15} '.format(item.name) +
                     '{0:<65}'.format(item.__repr__()))
        s.append(hr)
        return '\n'.join(s)


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
