import logging
import time
import uuid
import openpnm
import numpy as np
from copy import deepcopy
from openpnm.utils import HealthDict, Workspace
from openpnm.utils import SettingsAttr
from ._grid import Tableist


ws = Workspace()
logger = logging.getLogger(__name__)

__all__ = ['Project']


class ProjectSettings(SettingsAttr):
    r"""
    uuid : str
        A universally unique identifier for the object to keep things straight
    """
    uuid = ''


class Project(list):
    r"""
    This class provides a container for all OpenPNM objects in a given
    simulation.

    A simulation is defined as a Network and all of it's associated objects.
    When instantiating a Network, a Project can be passed as an argument, but
    if not given one is created. When instantiating any other object either
    a Network or a Project can be supplied. In the former case, the
    Network's Project is retrieved and used. The end result is that all
    objects are stored in a specific Project.

    The Project to which any object belongs can be retrieved with
    ``obj.project``. Conversely, printing a Project displays a list of all
    objects it contains.

    Moreover, all Projects are registered with the Workspace. Since there can
    be only instance of the Workspace it is possible to view all open Projects
    by printing the Workspace.

    See Also
    --------
    Workspace

    """

    def __init__(self, *args, **kwargs):
        name = kwargs.pop('name', None)
        super().__init__(*args, **kwargs)
        self.settings = ProjectSettings()
        ws[name] = self  # Register self with workspace
        self.settings['uuid'] = str(uuid.uuid4())

    def extend(self, obj):
        r"""
        This function is used to add objects to the project. Arguments can
        be single OpenPNM objects, an OpenPNM project list, or a plain list of
        OpenPNM objects. Note that if an object has the same name as one
        already existing on the project, the it will be renamed automatically.
        """
        if not isinstance(obj, list):
            obj = [obj]
        for item in obj:
            if hasattr(item, '_mro'):
                if 'GenericNetwork' in item._mro():
                    if self.network:
                        raise Exception('Project already has a network')
                # Must use append since extend breaks the dicts up into
                # separate objects, while append keeps it as a single object.
                if item in self:
                    raise Exception('Supplied object already part of project')
                if item.name in self.names:
                    item.name = self._generate_name(item)
                super().append(item)
            else:
                raise Exception('Only OpenPNM objects can be added')

    def append(self, obj):
        r"""
        The Project (a list) must be kept as a flat list, so the append
        function, which can normally be used to insert a list into a list, is
        overloaded to basically prevent the normal append operation and simply
        calls ``extend``.

        """
        self.extend(obj)

    def remove(self, obj):
        r"""
        The given object is removed from the project

        This removes the object, along with all it's labels in associated
        objects, but does NOT remove the associated objects.

        See Also
        --------
        purge_object

        """
        self.purge_object(obj, deep=False)

    def pop(self, index):
        r"""
        The object at the given index is removed from the list and returned.

        Notes
        -----
        This method uses ``purge_object`` to perform the actual removal of the
        object. It is reommended to just use that directly instead.

        See Also
        --------
        purge_object

        """
        obj = self[index]
        self.purge_object(obj, deep=False)
        return obj

    def insert(self, index, obj):
        r"""
        Inserts the supplied object at the specified index in the Project list

        Notes
        -----
        The order of the objects in an OpenPNM Project lists do not matter, so
        it is recommended to just use ``append`` instead.

        See Also
        --------
        append
        extend

        """
        self.extend(obj)

    def clear(self, objtype=[]):
        r"""
        Clears objects from the project entirely or selectively, depdening on
        the received arguments.

        Parameters
        ----------
        objtype : list[str]
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

    def copy(self, name=None):
        r"""
        Creates a deep copy of the current project

        A deep copy means that new, unique versions of all the objects are
        created but with identical data and properties.

        Parameters
        ----------
        name : str
            The name to give to the new project. If not supplied, a name
            is automatically generated.

        Returns
        -------
        proj : list
            A new Project object containing copies of all objects

        Notes
        -----
        Because they are new objects, they are given a new uuid
        (``obj.settings['_uuid']``), but the uuid of the original object
        is also stored (``obj.settings['_uuid_old']``) for reference.

        """
        if name is None:
            name = ws._gen_name()
        proj = deepcopy(self)
        for item in proj:
            item.settings['_uuid'] = str(uuid.uuid4())
        self.settings['_uuid'] = str(uuid.uuid4())
        ws[name] = proj
        return proj

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
        if isinstance(key, str):
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
        r"""
        Finds the Phase associated with a given object.

        Parameters
        ----------
        obj : Base
            Can either be a Physics or Algorithm object

        Returns
        -------
        phase : GenericPhase

        Raises
        ------
        If no Phase object can be found, then an Exception is raised.

        """
        # If received algorithm, tell user to use alg.settings.phase
        if obj._isa('algorithm'):
            raise Exception("The name(s) of the associated phase(s) with algorithm"
                            " objects can be found via 'alg.settings.phase'")
        # If received phase, just return self
        if obj._isa('phase'):
            return obj
        # Otherwise find it using bottom-up approach (i.e. look in phase keys)
        for item in self.phases().values():
            if (f'pore.{obj.name}' in item) or (f'throat.{obj.name}' in item):
                return item
        # If all else fails, throw an exception
        raise Exception(f'Cannot find a phase associated with {obj.name}')

    def find_geometry(self, physics):
        r"""
        Finds the Geometry associated with a given Physics

        Parameters
        ----------
        physics : GenericPhysics
            The Physics object for which the Geometry object is sought

        Returns
        -------
        geom : GenericGeometry

        Raises
        ------
        If no Geometry object can be found, then an Exception is raised.

        """
        # Otherwise, use the bottom-up approach
        for geo in self.geometries().values():
            if physics in self.find_physics(geometry=geo):
                return geo
        # If all else fails, throw an exception
        raise Exception('Cannot find a geometry associated with '+physics.name)

    def find_physics(self, geometry=None, phase=None):
        r"""
        Finds the Physics object(s) associated with a given Geometry,
        Phase, or combination.

        Parameters
        ----------
        geometry : GenericGeometry
            The Geometry object for which the Physics object(s) are sought
        phase : GenericPhase
            The Phase object for which the Physics object(s) are sought

        Returns
        -------
        physics : list
            A list containing the Physics object(s). If only a
            ``geometry`` is specified the the Physics for all Phases is
            returned. If only a ``phase`` is specified, then the Physics
            for all Geometries is returned.  If both ``geometry`` and
            ``phase`` is specified then the list only contains a single
            Physics. If no Physics is found, the the list will be empty.
            See the Notes section for more information.

        See Also
        --------
        grid

        Notes
        -----
        The Project has a ``grid`` attribute that shows the associations of
        all objects.  If each Geometry represents a row and each Phase is a
        column, then each row/col intersection represents a Physics. This
        method finds the Physics' at each intersection

        """

        if geometry is not None and phase is not None:
            physics = self.find_physics(geometry=geometry)
            phases = list(self.phases().values())
            return physics[phases.index(phase)]

        if geometry is not None and phase is None:
            result = []
            geoPs = self.network[f'pore.{geometry.name}']
            geoTs = self.network[f'throat.{geometry.name}']
            for _phase in self.phases().values():
                physics = self.find_physics(phase=_phase)
                for phys in physics:
                    Ps = _phase.pores(phys.name)
                    physPs = _phase.to_mask(pores=Ps)
                    Ts = _phase.throats(phys.name)
                    physTs = _phase.to_mask(throats=Ts)
                    if np.all(geoPs == physPs) and np.all(geoTs == physTs):
                        result.append(phys)
            return result

        if geometry is None and phase is not None:
            names = set(self.physics().keys())
            keys = set([item.split('.')[-1] for item in phase.keys()])
            hits = names.intersection(keys)
            return [self.physics().get(i, None) for i in hits]

        return list(self.physics().values())

    def find_full_domain(self, obj):
        r"""
        Finds the full domain object associated with a given object.
        For geometry the network is found, for physics the phase is found
        and for all other objects which are defined for for the full
        domain, themselves are found.

        Parameters
        ----------
        obj : Base
            Can be any object

        Returns
        -------
        obj : Base

        """
        if 'Subdomain' not in obj._mro():
            # Network, Phase, Algorithm
            return obj
        if obj._isa() == 'geometry':
            # Geometry
            return self.network
        # Physics
        return self.find_phase(obj)

    def _validate_name(self, name):
        if name in self.names:
            raise Exception('Another object is already named '+name)
        for item in self:
            for key in item.keys():
                if key.split('.')[1] == name:
                    raise Exception('A property/label is already named '+name)

    def _generate_name(self, obj):
        prefix = obj.settings['prefix']
        num = len(self._get_objects_by_type(obj._isa())) + 1
        name = prefix + '_' + str(num).zfill(2)
        try:
            self._validate_name(name)
        except Exception:
            name = prefix + '_' + str(np.random.randint(100, 999))
        return name

    @property
    def names(self):
        names = [i.name for i in self]
        return names

    def purge_object(self, obj, deep=False):
        r"""
        Removes an object from the Project.  This removes all references
        to the object from all other objects (i.e. removes labels)

        Parameters
        ----------
        obj : Base or list[Base]
            The object(s) to purge
        deep : bool
            A flag that indicates whether to remove associated objects.
            If ``True``, then removing a Geometry or Phase also removes
            the associated Physics objects. If ``False`` (default) then
            only the given object is removed, along with its labels in all
            associated objects. Removing a Physics always keeps associated
            Geometry and Phases since they might also be associated with
            other Physics objects.

        Returns
        -------
        None

        Raises
        ------
        An Exception is raised if the object is a Network.

        Notes
        -----
        For a clearer picture of this logic, type ``print(project.grid)`` at
        the console.  A deep purge of a Geometry is like removing a row, while
        a Phase is like removing a column.

        """
        if isinstance(obj, list):
            for item in obj:
                self.purge_object(obj=item, deep=deep)
            return
        if obj._isa() in ['physics', 'algorithm']:
            self._purge(obj)
        if obj._isa() == 'geometry':
            if deep:
                physics = self.find_physics(geometry=obj)
                for phys in physics:
                    self._purge(self.physics()[phys.name])
            self._purge(obj)
        if obj._isa() == 'phase':
            if deep:
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
        super().remove(obj)

    def save_object(self, obj):
        r"""
        Saves the given object or list of objects to a pickle file

        Parameters
        ----------
        obj : Base or list[Base]
            The objects to be saved. Depending on the object type, the
            file extension will be one of 'net', 'geo', 'phase', 'phys'
            or 'alg'.

        Returns
        -------
        None

        """
        from openpnm.io import Pickle
        Pickle.save_object_to_file(objs=obj)

    def load_object(self, filename):
        r"""
        Loads a single object from a pickle file

        Parameters
        ----------
        filename : str or Path
            The name of the file containing the saved object.  Can include
            an absolute or relative path as well. If only a filename is
            given it will be saved in the current working directory. The
            object type is inferred from

        Returns
        -------
        None

        """
        from openpnm.io import Pickle
        Pickle.load_object_from_file(filename=filename, project=self)

    def save_project(self, filename=None):
        r"""
        Saves the current project to a ``pnm`` file.

        Parameters
        ----------
        filename : str or Path
            The name of the file. Can include an absolute or relative path
            as well. If only a filename is given it will be saved in the
            current working directory.

        Returns
        -------
        None

        """
        ws.save_project(project=self, filename=filename)

    def _new_object(self, objtype, name=None):
        if objtype.startswith('net'):
            obj = openpnm.network.GenericNetwork(project=self, name=name)
        elif objtype.startswith('geo'):
            obj = openpnm.geometry.GenericGeometry(network=self.network, name=name,
                                                   pores=[], throats=[])
        elif objtype.startswith('pha'):
            obj = openpnm.phase.GenericPhase(network=self.network, name=name)
        elif objtype.startswith('phy'):
            obj = openpnm.physics.GenericPhysics(network=self.network, name=name)
        elif objtype.startswith('alg'):
            obj = openpnm.algorithms.GenericAlgorithm(network=self.network, name=name)
        else:
            obj = openpnm.core.Base(project=self, name=name)
        return obj

    def export_data(self, phases=[], filename=None, filetype=None):
        r"""
        Exports the pore and throat data from the given object(s) into
        the specified file and format.

        Parameters
        ----------
        phases : list[GenericPhase]
            The data on each supplied phase will be added to file
        filename : str
            The file name to use.  If none is supplied then one will be
            automatically generated based on the name of the project
            containing the supplied Network, with the date and time appended.
        filetype : str
            Which file format to store the data. If a valid extension is
            included in the ``filename``, this is ignored. Option are:

            ======== ==========================================================
            format   description
            ======== ==========================================================
            'vtk'    (default) The Visualization Toolkit format, used by
                     various softwares such as Paraview. This actually
                     produces a 'vtp' file. NOTE: This can be quite slow since
                     all the data is written to a simple text file. For large
                     data simulations consider 'xdmf'.

            'csv'    The comma-separated values format, which is easily
                     opened in any spreadsheet program. The column names
                     represent the property name, including the type and name
                     of the object to which they belonged, all separated by
                     the pipe character.

            'xmf'    The extensible data markup format, is a very efficient
                     format for large data sets. This actually results in the
                     creation of two files, the *xmf* file and an associated
                     *hdf* file. The *xmf* file contains instructions for
                     looking into the *hdf* file where the data is stored.
                     Paraview opens the *xmf* format natively and is very fast.

            'pkl'    A python pickle files which is a saved represenation of
                     the entire simulation.

            'mat'    Matlab 'mat-file', which can be opened in Matlab.
            ======== ==========================================================

        Notes
        -----
        This is a helper function for the actual functions in the ``io``
        module. For more control over the format of the output, and more
        information about the format refer to ``openpnm.io``.

        """
        from openpnm import io

        project = self
        network = self.network
        if filename is None:
            filename = project.name + '_' + time.strftime('%Y%b%d_%H%M%p')
        if filetype is None:
            if '.' in filename:
                filetype = filename.split('.')[-1]
            else:
                raise Exception('File type not given')
        # Convert file type to io class name
        temp = {"hdf": "hdf5",
                "hdf5": "hdf5",
                "xmf": "xdmf",
                "xdmf": "xdmf",
                "vtk": "vtk",
                "vtp": "vtk",
                "pkl": "pickle",
                "mat": "matlab",
                "csv": "csv"}
        if filetype not in list(temp.keys()):
            raise Exception(f'File type {filetype} is not supported')
        if filetype in temp.keys():
            filetype = temp[filetype]

        # Fetch correct io class, using case insensitive look-up
        fmt = getattr(io, 'to_' + filetype)
        fmt(network=network, phases=phases, filename=filename)

    @property
    def network(self):
        net = list(self._get_objects_by_type('network').values())
        net = net[0] if len(net) > 0 else None
        return net

    def geometries(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('geometry')

    def phases(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('phase')

    def physics(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('physics')

    def algorithms(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('algorithm')

    def _get_object_by_name(self, name):
        for item in self:
            if item.name == name:
                return item
        raise Exception('An object named ' + name + ' was not found')

    def _get_objects_by_type(self, objtype):
        return {item.name: item for item in self if item._isa(objtype)}

    def __str__(self):
        s = []
        hr = '―'*78
        s.append(hr)
        s.append(' {0:<15} '.format('Object Name')
                 + '{0:<65}'.format('Object ID'))
        s.append(hr)
        for item in self:
            s.append(' {0:<15} '.format(item.name)
                     + '{0:<65}'.format(item.__repr__()))
        s.append(hr)
        return '\n'.join(s)

    def __repr__(self):
        return self.__str__()

    def check_geometry_health(self):
        r"""
        Performs a check to find pores with overlapping or undefined
        Geometries.

        Returns
        -------
        A HealthDict

        """
        health = HealthDict()
        health['overlapping_pores'] = []
        health['undefined_pores'] = []
        health['overlapping_throats'] = []
        health['undefined_throats'] = []
        geoms = self.geometries().keys()
        if len(geoms) > 0:
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
        else:
            health['undefined_pores'] = self.network.Ps
            health['undefined_throats'] = self.network.Ts
        return health

    def check_physics_health(self, phase):
        r"""
        Performs a check to find pores which have overlapping or missing
        Physics.

        Parameters
        ----------
        phase : GenericPhase
            The Phase whose Physics should be checked

        Returns
        -------
        HealthDict

        """
        health = HealthDict()
        health['overlapping_pores'] = []
        health['undefined_pores'] = []
        health['overlapping_throats'] = []
        health['undefined_throats'] = []
        geoms = self.geometries().keys()
        if len(geoms) > 0:
            phys = self.find_physics(phase=phase)
            if len(phys) == 0:
                raise Exception(str(len(geoms))+' geometries were found, but'
                                + ' no physics')
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

    def check_data_health(self, obj):
        r"""
        Checks the health of pore and throat data arrays.

        Parameters
        ----------
        obj : Base
            A handle of the object to be checked

        Returns
        -------
        health : dict
            Returns a HealthDict object which is a basic dictionary with an
            added ``health`` attribute that is ``True`` is all entries in the
            dict are deemed healthy (empty lists), or ``False`` otherwise.

        """
        health = HealthDict()
        for item in obj.props():
            health[item] = []
            if obj[item].dtype == 'O':
                health[item] = 'No checks on object'
            elif np.sum(np.isnan(obj[item])) > 0:
                health[item] = 'Has NaNs'
            elif np.shape(obj[item])[0] != obj._count(item.split('.')[0]):
                health[item] = 'Wrong Length'
        return health

    def check_network_health(self):
        r"""
        This method checks the topological health of the network

        The following aspects are checked for:

            (1) Isolated pores
            (2) Disconnected clusters of pores
            (3) Duplicate throats
            (4) Headless throats
            (5) Bidirectional throats

        Returns
        -------
        health : dict
            A dictionary containing the offending pores or throat numbers
            under each named key.

        Notes
        -----
        It also returns a list of which pores and throats should be trimmed
        from the network to restore health.  This list is a suggestion only,
        and is based on keeping the largest cluster and trimming the others.

        - Does not yet check for duplicate pores
        - Does not yet suggest which throats to remove
        - This is just a 'check' and does not 'fix' the problems it finds

        """
        import openpnm.models.network as mods

        health = HealthDict()
        net = self.network

        # Check for headless throats
        headless = mods.headless_throats(net)
        health['headless_throats'] = np.where(headless)[0].tolist()

        # Check for throats that loop back onto the same pore
        looped = mods.looped_throats(net)
        health['looped_throats'] = np.where(looped)[0].tolist()

        # Check for individual isolated pores
        isolated = mods.isolated_pores(net)
        health['isolated_pores'] = np.where(isolated)[0].tolist()

        # Check for separated clusters of pores
        size = mods.cluster_size(net)
        mx = np.max(size)
        health['disconnected_pores'] = np.where(size < mx)[0].tolist()

        # Check for duplicate throats
        dupes = mods.duplicate_throats(net)
        health['duplicate_throats'] = np.where(dupes)[0].tolist()

        # Check for bidirectional throats
        bidir = mods.bidirectional_throats(net)
        health['bidirectional_throats'] = np.where(bidir)[0].tolist()

        return health

    def show_model_dependencies(self, prop, obj):
        """Brief explanation of 'show_model_dependencies'"""
        deps = {prop: self._get_deps(prop, obj)}
        self._view_dependencies(deps)

    def _get_deps(self, prop, obj):

        deps = {}
        try:
            model = obj.models[prop]
            for item in model.values():
                if isinstance(item, str):
                    if item.startswith('pore.') or item.startswith('throat.'):
                        upstream = self._get_deps(item, obj)
                        deps.update({item: upstream})
        except KeyError:
            if obj._isa('physics'):
                phase = self.find_phase(obj)
                geom = self.find_geometry(obj)
                if prop in phase.models.keys():
                    deps.update(self._get_deps(prop, phase))
                elif prop in geom.models.keys():
                    deps.update(self._get_deps(prop, geom))
                else:
                    pass
        return deps

    def _deps_to_jsongraph(self, children, name=None, parent=None):
        if parent is None:
            parent = "null"
        if name is None:
            name = list(children.keys())[0]
        tree = {"name": name,
                "parent": parent,
                "color": hex(hash(name.split('.')[1]))[3:9],
                "children": []}
        for item in children[name].keys():
            sub_tree = self._deps_to_jsongraph(parent=name, name=item,
                                               children=children[name])
            tree["children"].append(sub_tree)
        return tree

    def _view_dependencies(self, deps, port=8008):
        import json
        import webbrowser
        import threading
        import os
        web_dir = os.path.join(os.path.dirname(__file__), '../../public')
        os.chdir(web_dir)
        from http.server import HTTPServer, SimpleHTTPRequestHandler
        server = HTTPServer(server_address=('', port),
                            RequestHandlerClass=SimpleHTTPRequestHandler)
        thread = threading.Thread(target=server.serve_forever)
        thread.daemon = True
        thread.start()

        data = self._deps_to_jsongraph(deps)
        with open('data/tree.json', 'w') as outfile:
            json.dump(data, outfile)

        # Launch browser
        webbrowser.open(f"http://localhost:{port}/dep_map.html")

    def inspect_locations(self, element, indices, objs=[], mode='all'):
        r"""
        Shows the values of all props and/or labels for a given subset of
        pores or throats.

        Parameters
        ----------
        element : str
            The type of locations to inspect, either 'pores', or 'throats'
        indices : array_like
            The pore or throat indices to inspect
        objs : list[Base]
            If given, then only the properties on the recieved object are
            inspected.  If not given, then all objects are inspected (default).
        mode : list[str]
            Indicates whether to inspect 'props', 'labels', or 'all'.  The
            default is all

        Returns
        -------
        df : Pandas DataFrame
            A data frame object with each location as a column and each row
            as a property and/or label.

        """
        from pandas import DataFrame
        props = {}
        if not isinstance(objs, list):
            objs = [objs]
        if len(objs) == 0:
            objs = self
        for obj in objs:
            d = {k: obj[k][indices] for k in obj.keys(element=element, mode=mode)}
            for item in list(d.keys()):
                if d[item].ndim > 1:
                    d.pop(item)
                    if item == 'pore.coords':
                        d['pore.coords_X'], d['pore.coords_Y'], \
                            d['pore.coords_Z'] = obj['pore.coords'][indices].T
                    if item == 'throat.conns':
                        d['throat.conns_head'], d['throat.conns_tail'] = \
                            obj['throat.conns'][indices].T
                _ = [props.update({obj.name+'.'+item: d[item]}) for item in d.keys()]
        df = DataFrame(props)
        df = df.rename(index={k: indices[k] for k, _ in enumerate(indices)})
        return df.T

    def _regenerate_models(self, objs=[], propnames=[]):
        r"""
        Can be used to regenerate models across all objects in the project.

        Parameters
        ----------
        objs : list[Base]
            Can be used to specify which specific objects to regenerate. The
            default is to regenerate all objects. If a subset of objects is
            given, this function ensure they are generated in a sensible order
            such as any phases are done before any physics objects.
        propnames : list[str], or str
            The specific model to regenerate. If none are given then ALL
            models on all objects are regenerated. If a subset is given,
            then only object that have a corresponding model are regenerated,
            to avoid any problems. This means that a single model can be
            given, without specifying the objects.

        """
        objs = list(objs)
        if objs == []:
            objs = self
        if isinstance(propnames, str):
            propnames = [propnames]
        # Sort objs in the correct order (geom, phase, phys)
        net = [i for i in objs if i is self.network]
        geoms = [i for i in objs if i in self.geometries().values()]
        phases = [i for i in objs if i in self.phases().values()]
        phys = [i for i in objs if i in self.physics().values()]
        objs = net + geoms + phases + phys
        for obj in objs:
            if len(propnames) > 0:
                for model in propnames:
                    if model in obj.models.keys():
                        obj.regenerate_models(propnames=model)
            else:
                obj.regenerate_models()

    def _generate_grid(self):
        grid = ProjectGrid()
        # Create first/index column of grid
        rows = [self.network.name] + list(self.geometries().keys())
        grid.add_row(num=len(rows) - 1)
        grid.set_col(0, rows)
        # Creatle first/header row of grid
        cols = list(self.phases().keys())
        grid.add_col(num=len(cols))
        grid.set_row(0, vals=[self.network.name] + cols)
        # Now add physics objects to grid, adding new columns/rows as needed.
        miss = 0
        for p in self.physics().values():
            try:
                row = self.find_geometry(p)
                try:
                    col = self.find_phase(p)
                    grid.set_row_and_col(row=row.name, col=col.name, val=p.name)
                except Exception:
                    miss += 1
                    grid.set_row_and_col(row=row.name, col='?'*miss, val=p.name)
            except:
                try:
                    col = self.find_phase(p)
                    miss += 1
                    grid.set_row_and_col(row='?'*miss, col=col.name, val=p.name)
                except Exception:
                    miss += 1
                    grid.set_row_and_col(row='?'*miss, col='?'*miss, val=p.name)
        # See if any pores/throats are not assigned and add blank row
        if len(self.geometries()) > 0:
            h = self.check_geometry_health()
            if (len(h['undefined_pores']) > 0) or (len(h['undefined_throats']) > 0):
                grid.add_row()
        return grid

    def _get_grid(self):
        if not hasattr(self, '_grid'):
            grid = self._generate_grid()
            self._grid = grid
        else:  # Update current grid with new data, to save formats and settings
            grid = self._generate_grid()
            self._grid._grid.table_data = grid._grid.table_data
        return self._grid

    def _set_grid(self, grid):
        self._grid = grid

    grid = property(fget=_get_grid, fset=_set_grid)


class ProjectGrid(Tableist):
    r"""
    This is a subclass of a Tableist grid, which adds the ability to lookup
    by geometries and phases, as more specific versions of rows and cols
    """

    def row(self, name):
        r"""
        Retrieve a specified row from the table

        Parameters
        ----------
        name : str
            The row name, specified by the ``geometry`` object name

        Returns
        -------
        table
            A table object containing only a single row

        """
        return self.get_row(name)._grid.table_data[0]

    def col(self, name):
        r"""
        Retrieves a specified column from the table

        Parameters
        ----------
        name : str
            The column name, specified by the ``phase`` object name

        Returns
        -------
        table
            A table object containing only a single column

        """
        temp = self.get_col(name)._grid.table_data
        temp = [i[0] for i in temp]
        return temp

    def geometries(self):
        r"""
        Retrieves a list of all geometries
        """
        temp = self.index[1:]
        temp = [i[0] for i in temp]
        return temp

    def phases(self):
        r"""
        Retrieves a list of all phases
        """
        return self.header[0][1:]
