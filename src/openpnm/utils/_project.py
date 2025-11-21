import logging
import pickle
import re
import uuid
from copy import deepcopy
from datetime import datetime

import numpy as np

from openpnm.utils import SettingsAttr, Workspace

ws = Workspace()
logger = logging.getLogger(__name__)


__all__ = [
    'Project',
]


class ProjectSettings(SettingsAttr):
    r"""
    uuid : str
        A universally unique identifier for the object to keep things straight
    """
    uuid = ''
    original_uuid = ''

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        name = ws._validate_name(name)
        for v in list(ws.values()):
            if v.settings is self:
                ws[name] = ws.pop(v.settings.name)
        self._name = name


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
        self.settings = ProjectSettings()
        self.settings['name'] = kwargs.pop('name', None)
        self.settings['uuid'] = str(uuid.uuid4())
        self.settings['original_uuid'] = self.settings['uuid']
        super().__init__(*args, **kwargs)
        ws[self.settings['name']] = self

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except TypeError:
            if isinstance(key, str):  # Enable dict-style lookup if key is a string
                for item in self:
                    if item.name == key:
                        return item
            raise KeyError(key)

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
        (``obj.settings['uuid']``), but the uuid of the original object
        is also stored (``obj.settings['original_uuid']``) for reference.

        """
        name = ws._validate_name(name)
        proj = deepcopy(self)
        for item in proj:
            item.settings['uuid'] = str(uuid.uuid4())
        proj.settings['uuid'] = str(uuid.uuid4())
        proj.settings['name'] = name
        ws[name] = proj
        return proj

    def _set_name(self, name):
        self.settings['name'] = name

    def _get_name(self):
        return self.settings['name']

    name = property(fget=_get_name, fset=_set_name)

    def _validate_name(self, name):
        if name in self.names:
            raise Exception('Another object is already named '+name)
        for item in self:
            for key in item.keys():
                if key.split('.')[1] == name:
                    raise Exception('A property/label is already named '+name)

    def _generate_name(self, name=''):
        if name in [None, '']:
            name = 'obj_01'  # Give basic name, then let rest of func fix it
        warn = True
        if name.endswith('_?'):
            name = name.replace('_?', '_01')
            warn = False
        if name in self.names:  # If proposed name is taken, increment it
            proposed_name = name
            if not re.search(r'_\d+$', name):  # If name does not end with _##
                name = name + '_01'
            prefix, count = name.rsplit('_', 1)
            n = [0]
            for item in self:
                if item.name.startswith(prefix+'_'):
                    n.append(int(item.name.split(prefix+'_')[1]))
            name = prefix+'_'+str(max(n)+1).zfill(2)
            if warn:
                logger.warn(f'{proposed_name} is already taken, using {name} instead')
        self._validate_name(name)
        return name

    @property
    def names(self):
        names = [i.name for i in self]
        return names

    @property
    def network(self):
        for item in self:
            if ('throat.conns' in item.keys()) or ('pore.coords' in item.keys()):
                return item

    @property
    def phases(self):
        from openpnm.phase import Phase
        phases = []
        for item in self:
            if isinstance(item, Phase):
                phases.append(item)
        return phases

    @property
    def algorithms(self):
        from openpnm.algorithms import Algorithm
        algs = []
        for item in self:
            if isinstance(item, Algorithm):
                algs.append(item)
        return algs

    @property
    def workspace(self):
        return ws

    def _get_locations(self, label):
        r"""
        Find locations indicated by the given label regardless of which object
        it is defined on

        Parameters
        ----------
        label : str
            The label whose locations are sought, such as 'pore.left'

        Returns
        -------
        locations : ndarray
            A boolean array with ``True`` values indicating which locations
            have the given label

        Notes
        -----
        The returns the first instance of ``label`` that it finds
        """
        for item in self:
            try:
                return item[label]
            except KeyError:
                pass
        raise KeyError(label)

    def __str__(self):  # pragma: no cover
        hr = '―'*78
        s = '═'*78 + '\n'
        s += 'Object Name : Object Class and ID' + '\n'
        s += hr + '\n'
        for item in self:
            s += item.__repr__() + '\n'
        s += hr
        return s
