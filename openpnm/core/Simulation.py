import time
import pandas as pd
from openpnm.core import Workspace
ws = Workspace()


class Simulation(dict):

    def __init__(self, name=None):
        super().__init__()
        if name is None:
            name = 'sim'+str(len(ws.keys())+1).zfill(3)
        self.name = name
        self.grid = Grid()
        self.pore_map = pd.DataFrame()
        self.throat_map = pd.DataFrame()
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

    def add_network(self, network):
        if self.network is not None:
            raise Exception("This simulation already has a Network")
        self.update({network.name: network})

    def add_phase(self, phase):
        if phase in self.phases.values():
            raise Exception("A Phase with that name has already been added")
        self.update({phase.name: phase})
        for geo in self.geometries.values():
            self.grid[geo.name][phase.name] = ''

    def add_geometry(self, geometry):
        if geometry in self.geometries.values():
            raise Exception("A Geometry with that name has already been added")
        self.update({geometry.name: geometry})
        self.grid[geometry.name] = {phase: '' for phase in self.phases}
        self.pore_map[geometry.name] = self.network['pore.'+geometry.name]
        self.throat_map[geometry.name] = self.network['throat.'+geometry.name]

    def add_physics(self, physics, geometry, phase):
        if physics in self.physics.values():
            raise Exception("A Physics with that name has already been added")
        self.update({physics.name: physics})
        self.grid[geometry.name][phase.name] = physics.name
        self.pore_map[physics.name] = phase['pore.'+physics.name]
        self.throat_map[physics.name] = phase['throat.'+physics.name]

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


class Grid(dict):

    def _get_phases(self):
        for sim in ws.values():
            if sim.grid is self:
                return list(sim.phases.keys())
    phases = property(fget=_get_phases)

    def _get_net(self):
        for sim in ws.values():
            if sim.grid is self:
                return sim.network
    network = property(fget=_get_net)

    def __repr__(self):
        s = []
        hr = '―'*80
        s.append(hr)
        fmt = ["| {"+str(i)+":^13} " for i in range(len(self.phases))]
        phases = [item for item in self.phases]
        s.append('| {0:^13}'.format(self.network.name) +
                 ''.join(fmt).format(*phases) + '|')
        s.append(hr)
        for geo in self.keys():
            ind = '| {0:^13}'.format(geo)
            row = list(self[geo].values())
            s.append(ind + ''.join(fmt).format(*row) + '|')
            s.append(hr)
        return '\n'.join(s)
