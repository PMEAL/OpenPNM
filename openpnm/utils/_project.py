import logging
import uuid
import numpy as np
from copy import deepcopy
from openpnm.utils import Workspace
from openpnm.utils import SettingsAttr


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

    def __getitem__(self, key):
        if isinstance(key, str):  # Enable dict-style retrieval if key is a string
            obj = None
            for item in self:
                if item.name == key:
                    obj = item
            if obj is None:
                raise KeyError(key)
        else:
            obj = super().__getitem__(key)
        return obj

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
            # This needs to append or else it unpacks the dict into arrays
            super().append(item)

    def append(self, obj):
        self.extend(obj)

    def insert(self, index, obj):
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

    @property
    def network(self):
        for item in self:
            if 'throat.conns' in item.keys():
                return item

    def find_phase(self, obj):
        return obj

    def find_full_domain(self, obj):
        return obj

    def phases(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('phase')

    def algorithms(self, name=None):
        if name:
            return self._get_object_by_name(name)
        return self._get_objects_by_type('algorithm')

    def __str__(self):
        s = []
        hr = '―'*78
        s.append('═'*78)
        s.append(' {0:<15} '.format('Object Name')
                 + '{0:<65}'.format('Object ID'))
        s.append(hr)
        for item in self:
            s.append(' {0:<15} '.format(item.name)
                     + '{0:<65}'.format(item.__repr__()))
        s.append(hr)
        return '\n'.join(s)

    # def __repr__(self):
        # return self.__str__()
