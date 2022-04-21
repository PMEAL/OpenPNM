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
            super().append(item)

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

    def _validate_name(self, name):
        if name in self.names:
            raise Exception('Another object is already named '+name)
        for item in self:
            for key in item.keys():
                if key.split('.')[1] == name:
                    raise Exception('A property/label is already named '+name)

    def _generate_name(self, obj):
        prefix = obj.settings['prefix']
        num = len(self) + 1
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

    def purge_object(self, obj):
        r"""
        Removes an object from the Project.  This removes all references
        to the object from all other objects (i.e. removes labels)

        Parameters
        ----------
        obj : Base or list[Base]
            The object(s) to purge

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
                self.purge_object(obj=item)
            return
        for item in self:
            for key in list(item.keys()):
                if key.split('.')[-1] == obj.name:
                    del item[key]
        super().remove(obj)

    @property
    def network(self):
        for item in self:
            if 'throat.conns' in item.keys():
                return item

    def phases(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('phase')

    def algorithms(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('algorithm')

    def _get_object_by_name(self, name):
        for item in self:
            if item.name == name:
                return item
        raise Exception('An object named ' + name + ' was not found')

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
