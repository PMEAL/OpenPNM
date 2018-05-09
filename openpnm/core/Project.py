import time
import pickle
import h5py
from openpnm.core import Workspace
from openpnm.utils.misc import SettingsDict, HealthDict, PrintableList
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
        self.comments = 'Using OpenPNM ' + openpnm.__version__

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

    def names(self):
        r"""
        Returns a list of all object names
        """
        names = PrintableList([obj.name for obj in self])
        return names

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
            if physics in self.find_physics(geometry=geo):
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
            phys = list(self.physics().values())
            return phys

    def _validate_name(self, name):
        if name in self.names:
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

    @property
    def names(self):
        names = [i.name for i in self]
        return names

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

    def import_data(self, filename):
        r"""
        """
        raise NotImplementedError('Use the io module to import data')

    def export_data(self, network=None, phases=[], filename=None,
                    filetype='vtp'):
        r"""
        Export the pore and throat data from the given object(s) into the
        specified file and format.

        Parameters
        ----------
        network: OpenPNM Network Object
            The network containing the data to be stored

        phases : list of OpenPNM Phase Objects
            The data on each supplied phase will be added to file

        filename : string
            The file name to use.  If none is supplied then one will be
            automatically generated based on the name of the project
            containing the supplied Network, with the date and time appended.

        filetype : string
            Which file format to store the data.  If a valid extension is
            included in the ``filename``, this is ignored.  Option are:

            **'vtk'** : (default) The Visualization Toolkit format, used by
            various softwares such as Paraview.  This actually produces a 'vtp'
            file.  NOTE: This can be quite slow since all the data is written
            to a simple text file.  For large data simulations consider
            'xdmf'.

            **'csv'** : The comma-separated values format, which is easily
            openned in any spreadsheet program.  The column names represent
            the property name, including the type and name of the object to
            which they belonged, all separated by the pipe character.

            **'xmf'** : The extensible data markup format, is a very efficient
            format for large data sets.  This actually results in the creation
            of two files, the *xmf* file and an associated *hdf* file.  The
            *xmf* file contains instructions for looking into the *hdf* file
            where the data is stored. Paraview opens the *xmf* format natively,
            and is very fast.

            **'mat'** : Matlab 'mat-file', which can be openned in Matlab.

        Notes
        -----
        This is a helper function for the actual functions in the IO module.
        For more control over the format of the output, and more information
        about the format refer to ``openpnm.io``.

        """
        project = network.project
        if filename is None:
            filename = project.name + '_' + time.strftime('%Y%b%d_%H%M%p')
        # Infer filetype from extension on file name...if given
        if '.' in filename:
            exts = ['vtk', 'vtp', 'vtu', 'csv', 'xmf', 'xdmf', 'hdf', 'hdf5',
                    'h5', 'mat']
            if filename.split('.')[-1] in exts:
                filename, filetype = filename.rsplit('.', 1)
        if filetype.lower() in ['vtk', 'vtp', 'vtu']:
            openpnm.io.VTK.save(network=network, phases=phases,
                                filename=filename)
        if filetype.lower() == 'csv':
            openpnm.io.CSV.save(network=network, phases=phases,
                                filename=filename)
        if filetype.lower() in ['xmf', 'xdmf']:
            openpnm.io.XDMF.save(network=network, phases=phases,
                                 filename=filename)
        if filetype.lower() in ['hdf5', 'hdf', 'h5']:
            f = openpnm.io.HDF5.to_hdf5(network=network, phases=phases,
                                        filename=filename)
            f.close()
        if filetype.lower() == 'mat':
            openpnm.io.MAT.save(network=network, phases=phases,
                                filename=filename)

    def _dump_data(self, mode=['props']):
        r"""
        Dump data from all objects in project to an HDF5 file

        Parameters
        ----------
        mode : string or list of strings
            The type of data to be dumped to the HDF5 file.  Options are:

            **'props'** : Numerical data such as 'pore.diameter'.  The default
            is only 'props'.

            **'labels'** : Boolean data that are used as labels.  Since this
            is boolean data it does not consume large amounts of memory and
            probably does not need to be dumped.

        """
        with h5py.File(self.name + '.hdf5') as f:
            for obj in self:
                for key in list(obj.keys()):
                    tempname = obj.name + '|' + '_'.join(key.split('.'))
                    arr = obj[key]
                    if 'U' in str(obj[key][0].dtype):
                        pass
                    elif 'all' in key.split('.'):
                        pass
                    else:
                        f.create_dataset(name='/'+tempname, shape=arr.shape,
                                         dtype=arr.dtype, data=arr)
            for obj in self:
                obj.clear(mode=mode)

    def _fetch_data(self):
        r"""
        Retrieve data from an HDF5 file and place onto correct objects in the
        project
        """
        with h5py.File(self.name + '.hdf5') as f:
            # Reload data into project
            for item in f.keys():
                obj_name, propname = item.split('|')
                propname = propname.split('_')
                propname = propname[0] + '.' + '_'.join(propname[1:])
                self[obj_name][propname] = f[item]

    @property
    def network(self):
        net = list(self._get_objects_by_type('network').values())
        if len(net) > 0:
            net = net[0]
        else:
            net = None
        return net

    def geometries(self, name=None):
        if name:
            return self._get_object_by_name(name)
        else:
            return self._get_objects_by_type('geometry')

    def phases(self, name=None):
        if name:
            return self._get_object_by_name(name)
        else:
            return self._get_objects_by_type('phase')

    def physics(self, name=None):
        if name:
            return self._get_object_by_name(name)
        else:
            return self._get_objects_by_type('physics')

    def algorithms(self, name=None):
        if name:
            return self._get_object_by_name(name)
        else:
            return self._get_objects_by_type('algorithm')

    def _get_object_by_name(self, name):
        for item in self:
            if item.name == name:
                return item
        raise Exception('An object named ' + name + ' was not found')

    def _get_objects_by_type(self, objtype):
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
        health = HealthDict()
        health['overlapping_pores'] = []
        health['undefined_pores'] = []
        health['overlapping_throats'] = []
        health['undefined_throats'] = []
        geoms = self.geometries().keys()
        if len(geoms):
            net = self.network
            Ptemp = np.zeros((net.Np,))
            Ttemp = np.zeros((net.Nt,))
            for item in geoms:
                Pind = net['pore.'+item]
                Tind = net['throat.'+item]
                Ptemp[Pind] = Ptemp[Pind] + 1
                Ttemp[Tind] = Ttemp[Tind] + 1
            health['overlapping_pores'] = np.where(Ptemp > 1)[0].tolist()
            health['undefined_pores'] = np.where(Ptemp == 0)[0].tolist()
            health['overlapping_throats'] = np.where(Ttemp > 1)[0].tolist()
            health['undefined_throats'] = np.where(Ttemp == 0)[0].tolist()
        return health

    def check_physics_health(self, phase):
        r"""
        Perform a check to find pores which have overlapping or missing Physics
        """
        health = HealthDict()
        health['overlapping_pores'] = []
        health['undefined_pores'] = []
        health['overlapping_throats'] = []
        health['undefined_throats'] = []
        geoms = self.geometries().keys()
        if len(geoms):
            phys = self.find_physics(phase=phase)
            if len(phys) == 0:
                raise Exception(str(len(geoms))+' geometries were found, but' +
                                ' no physics')
            if None in phys:
                raise Exception('Undefined physics found, check the grid')
            Ptemp = np.zeros((phase.Np,))
            Ttemp = np.zeros((phase.Nt,))
            for item in phys:
                    Pind = phase['pore.'+item.name]
                    Tind = phase['throat.'+item.name]
                    Ptemp[Pind] = Ptemp[Pind] + 1
                    Ttemp[Tind] = Ttemp[Tind] + 1
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
