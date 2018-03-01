import time
from openpnm.core import Workspace
from openpnm.utils.misc import SettingsDict
import numpy as np
ws = Workspace()


class Project(list):

    def __init__(self, name=None):
        super().__init__()
        # Register self with workspace
        ws[name] = self
        self._grid = {}
        self.settings = SettingsDict()

    def append(self, obj):
        if hasattr(obj, '_mro'):
            if 'GenericNetwork' in obj._mro():
                if self.network:
                    raise Exception('Project already has a network')
            super().append(obj)
        else:
            raise Exception('Only OpenPNM objects can be added')

    @property
    def workspace(self):
        return ws

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
        mro = obj._mro()
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
            row = self.grid.row(geometry.name)
            phys = [self.physics().get(i, None) for i in row]
        elif phase:
            col = self.grid.col(phase.name)
            phys = [self.physics().get(i, None) for i in col]
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
        num = str(len([item for item in self if item._isa() == obj._isa()]))
        name = prefix + '_' + num.zfill(2)
        return name

    def purge_object(self, obj):
        r"""
        """
        if obj._isa() in ['geometry', 'physics', 'algorithm']:
            self._purge(obj)
        if obj._isa() == 'phase':
            physics = self.find_physics(phase=obj)
            for phys in physics:
                self._purge(self.physics()[phys.name])
            self._purge(obj)
        if obj._isa() == 'network':
            raise Exception('Cannot purge a network, just make a new project')

    def _purge(self, obj):
        for item in self:
            for key in list(item.keys()):
                if key.split('.')[-1] == obj.name:
                    del item[key]
        self.remove(obj)

    def _get_net(self):
        for item in self:
            if 'GenericNetwork' in item._mro():
                return item

    network = property(fget=_get_net)

    def geometries(self):
        _dict = {}
        for item in self:
            if 'GenericGeometry' in item._mro():
                _dict.update({item.name: item})
        return _dict

    def phases(self):
        _dict = {}
        for item in self:
            if 'GenericPhase' in item._mro():
                _dict.update({item.name: item})
        return _dict

    def physics(self):
        _dict = {}
        for item in self:
            if 'GenericPhysics' in item._mro():
                _dict.update({item.name: item})
        return _dict

    def algorithms(self):
        _dict = {}
        for item in self:
            if 'GenericAlgorithm' in item._mro():
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
        grid = {}
        for geo in self.geometries().keys():
            grid[geo] = {}
            for phase in self.phases().values():
                grid[geo][phase.name] = ''
                for phys in self.physics().keys():
                    if phys in [n.split('.')[1] for n in phase.keys()]:
                        if np.sum(net['pore.'+geo][phase.pores(phys)]) > 0:
                            grid[geo][phase.name] = phys
        grid = ProjectGrid(self.name, grid)
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

    def __init__(self, name='', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    def index(self):
        return list(self.keys())

    def header(self):
        d = []
        for item in self.keys():
            d.extend([i for i in self[item].keys()])
        return list(set(d))

    def row(self, name):
        return list(self[name].values())

    def col(self, name):
        col = []
        for row in self.index():
            col.append(self[row][name])
        return col

    def __str__(self):
        s = []
        hr = '―'*(16*(len(self.header())+1))
        s.append(hr)
        fmt = ["| {"+str(i)+":^13} " for i in range(len(self.header()))]
        cols = [item for item in self.header()]
        s.append('| {0:^13}'.format(self.name) +
                 ''.join(fmt).format(*cols) + '|')
        s.append(hr)
        for row in self.index():
            ind = '| {0:^13}'.format(row)
            temp = list(self[row].values())
            s.append(ind + ''.join(fmt).format(*temp) + '|')
            s.append(hr)
        return '\n'.join(s)


class ProjectGrid(Grid):

    def geometries(self):
        return self.index()

    def phases(self):
        return self.header()
