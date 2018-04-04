import time
import pickle
import h5py
from openpnm.core import Workspace
from openpnm.utils.misc import SettingsDict, HealthDict
import openpnm
import numpy as np
ws = Workspace()


class Project(list):

    def __init__(self, name=None):
        super().__init__()
        # Register self with workspace
        ws[name] = self
        self._grid = {}
        self.settings = SettingsDict()

    def extend(self, obj):
        if hasattr(obj, '_mro'):
            if 'GenericNetwork' in obj._mro():
                if self.network:
                    raise Exception('Project already has a network')
            super().append(obj)
        else:
            raise Exception('Only OpenPNM objects can be added')

    def append(self, obj):
        self.extend(obj)

    def clear(self, objtype=[]):
        r"""
        Clears objects from the project entirely or selectively, depdening on
        the received arguments.

        Parameters
        ----------
        objtype : list of strings
            A list containing the object type(s) to be removed.  If no types
            are specified, then all objects are removed.  To clear only objects
            of a specific type, use *'network'*, *'geometry'*, *'phase'*,
            *'physics'*, or *'algorithm'*.  It's also possible to use
            abbreviations, like *'geom'*.

        """
        if len(objtype) == 0:
            super().clear()
        else:
            names = [obj.name for obj in self]
            for name in names:
                try:
                    obj = self[name]
                    for t in objtype:
                        if obj._isa(t):
                            self.purge_object(obj)
                except KeyError:
                    pass

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
            obj = None
            for item in self:
                if item.name == key:
                    obj = item
            if obj is None:
                raise KeyError(key)
        else:
            obj = super().__getitem__(key)
        return obj

    def find_phase(self, obj):
        # If received phase, just return self
        if obj._isa('phase'):
            return obj
        # If phase happens to be in settings (i.e. algorithm), look it up
        if 'phase' in obj.settings.keys():
            phase = self.phases()[obj.settings['phase']]
            return phase
        # Otherwise find it using bottom-up approach (i.e. look in phase keys)
        for phase in self.phases().values():
            if ('pore.'+obj.name in phase) or ('throat.'+obj.name in phase):
                return phase
        # If all else fails, throw an exception
        raise Exception('Cannot find a phase associated with '+obj.name)

    def find_geometry(self, physics):
        # If geometry happens to be in settings, look it up directly
        if 'geometry' in physics.settings.keys():
            geom = self.geometries()[physics.settings['geometry']]
            return geom
        # Otherwise, use the bottom-up approach
        for geo in self.geometries().values():
            phys = self.find_physics(geometry=geo)
            if phys is physics:
                return geo
        # If all else fails, throw an exception
        raise Exception('Cannot find a geometry associated with '+physics.name)

    def find_physics(self, geometry=None, phase=None):
        if geometry and phase:
            physics = self.find_physics(geometry=geometry)
            phases = list(self.phases().values())
            phys = physics[phases.index(phase)]
            return phys
        elif geometry:
            result = []
            net = self.network
            geoPs = net['pore.'+geometry.name]
            geoTs = net['throat.'+geometry.name]
            for phase in self.phases().values():
                physics = self.find_physics(phase=phase)
                temp = None
                for phys in physics:
                    Ps = phase.map_pores(phys['pore._id'])
                    physPs = phase.tomask(pores=Ps)
                    Ts = phase.map_throats(phys['throat._id'])
                    physTs = phase.tomask(throats=Ts)
                    if np.all(geoPs == physPs) and np.all(geoTs == physTs):
                        temp = phys
                result.append(temp)
            return result
        elif phase:
            names = set(self.physics().keys())
            keys = set([item.split('.')[-1] for item in phase.keys()])
            hits = names.intersection(keys)
            phys = [self.physics().get(i, None) for i in hits]
            return phys
        else:
            phys = []
            for geom in self.geometries().values():
                phys.append(self.find_physics(geometry=geom))
            return phys

    def _validate_name(self, name):
        names = [i.name for i in self]
        if name in names:
            raise Exception('An object already exists named '+name)
        for item in self:
            for key in item.keys():
                if key.split('.')[1] == name:
                    raise Exception('A property/label is already named '+name)

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

    def save_object(self, obj):
        r"""
        """
        if not isinstance(obj, list):
            obj = [obj]
        for item in obj:
            filename = item.name + '.' + item.settings['prefix']
            pickle.dump({item.name: item}, open(filename, 'wb'))

    def load_object(self, filename):
        d = pickle.load(open(filename, 'rb'))
        for item in d.keys():
            self.extend(d[item])

    def _new_object(self, objtype, name=None):
        if objtype == 'network':
            obj = openpnm.network.GenericNetwork(project=self, name=name)
        elif objtype == 'geometry':
            obj = openpnm.geometry.GenericGeometry(project=self, name=name)
        elif objtype == 'phase':
            obj = openpnm.phases.GenericPhase(project=self, name=name)
        elif objtype == 'physics':
            obj = openpnm.physics.GenericPhysics(project=self, name=name)
        elif objtype == 'algorithm':
            obj = openpnm.algorithm.GenericAlgorithm(project=self, name=name)
        else:
            obj = openpnm.core.Base(project=self, name=name)
        return obj

    def dump_data(self):
        r"""
        Dump data from all objects in project to an HDF5 file
        """
        f = h5py.File(self.name + '.hdf5')
        try:
            for obj in self:
                for key in list(obj.keys()):
                    tempname = obj.name + '|' + '_'.join(key.split('.'))
                    if 'U' in str(obj[key][0].dtype):
                        pass
                    elif 'all' in key.split('.'):
                        pass
                    else:
                        arr = obj.pop(key)
                        f.create_dataset(name='/'+tempname, shape=arr.shape,
                                         dtype=arr.dtype, data=arr)
        except AttributeError:
            print('File is not empty, change project name and try again')
            f.close()
        f.close()

    def load_data(self):
        r"""
        Retrieve data from an HDF5 file and place onto correct objects in the
        project
        """
        f = h5py.File(self.name + '.hdf5')
        # Reload data into project
        for item in f.keys():
            obj_name, propname = item.split('|')
            propname = propname.split('_')
            propname = propname[0] + '.' + '_'.join(propname[1:])
            self[obj_name][propname] = f[item]
        f.close()

    @property
    def network(self):
        net = list(self._get_objects_of_type('network').values())
        if len(net) > 0:
            net = net[0]
        else:
            net = None
        return net

    def geometries(self):
        return self._get_objects_of_type('geometry')

    def phases(self):
        return self._get_objects_of_type('phase')

    def physics(self):
        return self._get_objects_of_type('physics')

    def algorithms(self):
        return self._get_objects_of_type('algorithm')

    def _get_objects_of_type(self, objtype):
        return {item.name: item for item in self if item._isa(objtype)}

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
        grid = {}
        row = {phase: '' for phase in self.phases().keys()}
        for geo in self.geometries().values():
            grid[geo.name] = row.copy()
            for phase in self.phases().values():
                phys = self.find_physics(phase=phase, geometry=geo)
                if phys is None:
                    phys = ''
                else:
                    phys = phys.name
                grid[geo.name][phase.name] = phys
        grid = ProjectGrid(self.network.name, grid)
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

    def check_geometry_health(self):
        r"""
        Perform a check to find pores with overlapping or undefined Geometries
        """
        geoms = self.geometries().keys()
        net = self.network
        Ptemp = np.zeros((net.Np,))
        Ttemp = np.zeros((net.Nt,))
        for item in geoms:
            Pind = net['pore.'+item]
            Tind = net['throat.'+item]
            Ptemp[Pind] = Ptemp[Pind] + 1
            Ttemp[Tind] = Ttemp[Tind] + 1
        health = HealthDict()
        health['overlapping_pores'] = np.where(Ptemp > 1)[0].tolist()
        health['undefined_pores'] = np.where(Ptemp == 0)[0].tolist()
        health['overlapping_throats'] = np.where(Ttemp > 1)[0].tolist()
        health['undefined_throats'] = np.where(Ttemp == 0)[0].tolist()
        return health

    def check_physics_health(self, phase):
        r"""
        Perform a check to find pores which have overlapping or missing Physics
        """
        phys = self.find_physics(phase=phase)
        if None in phys:
            raise Exception('Undefined physics found, check the grid')
        Ptemp = np.zeros((phase.Np,))
        Ttemp = np.zeros((phase.Nt,))
        for item in phys:
                Pind = phase['pore.'+item.name]
                Tind = phase['throat.'+item.name]
                Ptemp[Pind] = Ptemp[Pind] + 1
                Ttemp[Tind] = Ttemp[Tind] + 1
        health = HealthDict()
        health['overlapping_pores'] = np.where(Ptemp > 1)[0].tolist()
        health['undefined_pores'] = np.where(Ptemp == 0)[0].tolist()
        health['overlapping_throats'] = np.where(Ttemp > 1)[0].tolist()
        health['undefined_throats'] = np.where(Ttemp == 0)[0].tolist()
        return health


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
            s.append(ind + ''.join(fmt).format(*list(self[row].values()))+'|')
            s.append(hr)
        return '\n'.join(s)


class ProjectGrid(Grid):
    r"""
    This is a subclass of Grid, which adds the ability to lookup by geometries
    and phases, as more specific versions of rows and cols
    """

    def geometries(self):
        return self.index()

    def phases(self):
        return self.header()
