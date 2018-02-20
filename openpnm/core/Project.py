import time
from openpnm.core import Workspace
import numpy as np
ws = Workspace()


class Project(list):

    def __init__(self, name=None):
        super().__init__()
        # Register self with workspace
        ws[name] = self
        self._grid = {}
        self.settings = ws.settings.copy()

    def append(self, obj):
        if 'openpnm' in str(type(obj)):
            if 'GenericNetwork' in obj.mro():
                if self.network:
                    raise Exception('Project already has a network')
            super().append(obj)
        else:
            raise Exception('Only OpenPNM objects can be added')

    def _set_name(self, name):
        if name is None:
            name = ws._gen_name()
        ws[name] = self

    def _get_name(self):
        for key in ws.keys():
            if ws[key] is self:
                return key

    name = property(fget=_get_name, fset=_set_name)

    def __getitem__(self, key):
        if type(key) == str:
            for obj in self:
                if obj.name == key:
                    return obj
        else:
            return super().__getitem__(key)

    def find_phase(self, obj):
        mro = [c.__name__ for c in obj.__class__.__mro__]
        if 'GenericPhase'in mro:
            return obj
        if 'GenericAlgorithm' in mro:
            phase = self.phases()[obj.settings['phase']]
            return phase
        for g in self.geometries().values():
            for p in self.phases().values():
                if obj.name == self.grid[g.name][p.name]:
                    return p

    def find_geometry(self, physics):
        for g in self.geometries().values():
            for p in self.phases().values():
                if physics.name == self.grid[g.name][p.name]:
                    return g

    def find_physics(self, geometry=None, phase=None):
        if geometry and phase:
            name = self.grid[geometry.name][phase.name]
            phys = self[name]
        elif geometry:
            phys = self.grid.row(geometry)
        elif phase:
            phys = self.grid.col(phase)
        else:
            raise Exception('Must specify at least one of geometry or phase')
        if phys == ['']:
            phys = []
        return phys

    def _validate_name(self, name):
        names = [i.name for i in self]
        if name in names:
            raise Exception('An object with that name already exists!')
        for item in self:
            for key in item.keys():
                if key.split('.')[1] == name:
                    raise Exception('A property/label already uses that name')

    def _generate_name(self, obj):
        prefix = obj.settings['prefix']
        num = str(len([item for item in self if item.isa() == obj.isa()]))
#        if 'GenericNetwork' in obj.mro():
#            num = str(1).zfill(2)
#        elif 'GenericGeometry' in obj.mro():
#            num = str(len(self.geometries().keys())).zfill(2)
#        elif 'GenericPhase' in obj.mro():
#            num = str(len(self.phases().keys())).zfill(2)
#        elif 'GenericPhysics' in obj.mro():
#            num = str(len(self.physics().keys())).zfill(2)
#        elif 'GenericAlgorithm' in obj.mro():
#            num = str(len(self.algorithms().keys())).zfill(2)
#        else:
#            num = str(len(self)).zfill(2)
        name = prefix + '_' + num
        return name

    def _get_net(self):
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericNetwork' in mro:
                return item

    network = property(fget=_get_net)

    def geometries(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericGeometry' in mro:
                _dict.update({item.name: item})
        return _dict

    def phases(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericPhase' in mro:
                _dict.update({item.name: item})
        return _dict

    def physics(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericPhysics' in mro:
                _dict.update({item.name: item})
        return _dict

    def algorithms(self):
        _dict = {}
        for item in self:
            mro = [c.__name__ for c in item.__class__.__mro__]
            if 'GenericAlgorithm' in mro:
                _dict.update({item.name: item})
        return _dict

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
        for geo in self.geometries().keys():
            grid[geo] = {}
            for phase in self.phases().values():
                grid[geo][phase.name] = ''
                for phys in self.physics().keys():
                    if phys in [n.split('.')[1] for n in phase.keys()]:
                        if np.sum(net['pore.'+geo][phase.pores(phys)]) > 0:
                            grid[geo][phase.name] = phys
        self._grid = grid
        return grid

    grid = property(fget=_get_grid)

    def __str__(self):
        s = []
        hr = '―'*78
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
        return list(sim.geometries().keys())

    geometries = property(fget=_get_geometries)

    def _get_phases(self):
        sim = self._get_sim()
        return list(sim.phases().keys())

    phases = property(fget=_get_phases)

    def _get_net(self):
        sim = self._get_sim()
        return sim.network

    network = property(fget=_get_net)

    def row(self, geometry):
        return list(self[geometry.name].values())

    def col(self, phase):
        col = []
        for geo in self.geometries:
            col.append(self[geo][phase.name])
        return col

    def __str__(self):
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
