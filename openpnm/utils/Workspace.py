import pickle
import openpnm
import time
import copy
import warnings
import numpy as np
from pathlib import Path
from openpnm.utils import SettingsDict, logging
logger = logging.getLogger()


class SettingsDict(SettingsDict):
    def __setitem__(self, key, value):
        if key == 'loglevel':
            logger.setLevel(value)
        super().__setitem__(key, value)


class Workspace(dict):
    r"""
    The Workspace object provides the highest level of adminstrative
    control over active OpenPNM sessions.

    It is a ``dictionary`` that stores a list of all open Projects by name.

    This class is a
    `singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ so that
    whenever and where ever a Workspace is instantiated, the same instance is
    obtained.  This allows it to maintain a definitive record of all open
    Projects.

    See Also
    --------
    Project

    Notes
    -----
    The Workspace object contains a variety of functions that one might expect
    from the 'file-menu' in a typical GUI.  The table below provides a list
    along with a brief description.

    +----+------------------+-------------------------------------------------+
    | #  | Method           | Description                                     |
    +====+==================+=================================================+
    | 1  | new_project      | Creates a new empty Project object              |
    +----+------------------+-------------------------------------------------+
    | 2  | load_project     | Loads a Project from the specified 'pnm' file   |
    +----+------------------+-------------------------------------------------+
    | 3  | copy_project     | Make a copy of an existing Project              |
    +----+------------------+-------------------------------------------------+
    | 4  | load_workspace   | Loads all saved Projects in a 'pnm' file int... |
    +----+------------------+-------------------------------------------------+
    | 5  | save_project     | Saves given Project to a 'pnm' file             |
    +----+------------------+-------------------------------------------------+
    | 6  | save_workspace   | Saves all the current Projects to a 'pnm' file  |
    +----+------------------+-------------------------------------------------+
    | 7  | close_project    | Removes the specified Project from the Works... |
    +----+------------------+-------------------------------------------------+

    """

    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Workspace.__instance__ is None:
            Workspace.__instance__ = dict.__new__(cls)
        return Workspace.__instance__

    def __init__(self):
        super().__init__()
        self.settings = SettingsDict()

    def __setitem__(self, name, project):
        if name is None:
            name = self._gen_name()
        if name in self.keys():
            raise Exception("A project named " + name + " already exists")
        if project in self.values():
            self.pop(project.name, None)
        if not isinstance(project, openpnm.utils.Project):
            project = openpnm.utils.Project(project, name=name)
        super().__setitem__(name, project)

    def _create_console_handles(self, project):
        r"""
        Adds all objects in the given project to the console as variables
        with handle names taken from each object's name.
        """
        import __main__
        for item in project:
            __main__.__dict__[item.name] = item

    def save_workspace(self, filename=''):
        r"""
        Saves all the current Projects to a 'pnm' file

        Parameters
        ----------
        filename : string, optional
            If no filename is given, a name is genrated using the current
            time and date. See Notes for more information on valid file names.

        See Also
        --------
        save_project

        Notes
        -----
        The filename can be a string such as 'saved_file.pnm'.  The string can
        include absolute path such as 'C:\networks\saved_file.pnm', or can
        be a relative path such as '..\..\saved_file.pnm', which will look
        2 directories above the current working directory.  It can also be a
        path object object such as that produced by ``pathlib`` or
        ``os.path`` in the Python standard library.

        """
        if filename == '':
            filename = 'workspace' + '_' + time.strftime('%Y%b%d_%H%M%p')
        filename = self._parse_filename(filename=filename, ext='pnm')
        d = {}
        for sim in self.values():
            d[sim.name] = sim
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    def load_workspace(self, filename):
        r"""
        Loads all saved Projects in a 'pnm' file into the Workspace

        Any Projects present in the current Workspace will be deleted.

        Parameters
        ----------
        filename : string, optional
            The name of the file to open.  See Notes for more information.

        See Also
        --------
        load_project

        Notes
        -----
        The filename can be a string such as 'saved_file.pnm'.  The string can
        include absolute path such as 'C:\networks\saved_file.pnm', or can
        be a relative path such as '..\..\saved_file.pnm', which will look
        2 directories above the current working directory.  Can also be a
        path object object such as that produced by ``pathlib`` or
        ``os.path`` in the Python standard library.

        """
        self.clear()
        self.load_project(filename=filename)

    def save_project(self, project, filename=''):
        r"""
        Saves given Project to a 'pnm' file

        This will include all of associated objects, including algorithms.

        Parameters
        ----------
        project : OpenPNM Project
            The project to save.

        filename : string, optional
            If no filename is given, the given project name is used. See Notes
            for more information.

        See Also
        --------
        save_workspace

        Notes
        -----
        The filename can be a string such as 'saved_file.pnm'.  The string can
        include absolute path such as 'C:\networks\saved_file.pnm', or can
        be a relative path such as '..\..\saved_file.pnm', which will look
        2 directories above the current working directory.  Can also be a
        path object object such as that produced by ``pathlib`` or
        ``os.path`` in the Python standard library.

        """
        if filename == '':
            filename = project.name
        filename = self._parse_filename(filename=filename, ext='pnm')

        # Save dictionary as pickle
        d = {project.name: project}
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    def load_project(self, filename, overwrite=False):
        r"""
        Loads a Project from the specified 'pnm' file

        The loaded project is added to the Workspace . This will *not* delete
        any existing Projects in the Workspace and will rename any Projects
        being loaded if necessary.

        Parameters
        ----------
        filename : string or path object
            The name of the file to open.  See Notes for more information.

        See Also
        --------
        load_workspace

        Notes
        -----
        The filename can be a string such as 'saved_file.pnm'.  The string can
        include absolute path such as 'C:\networks\saved_file.pnm', or can
        be a relative path such as '..\..\saved_file.pnm', which will look
        2 directories above the current working directory.  Can also be a
        path object object such as that produced by ``pathlib`` or
        ``os.path`` in the Python standard library.

        """
        filename = self._parse_filename(filename=filename, ext='pnm')
        temp = {}  # Read file into temporary dict
        with open(filename, 'rb') as f:
            d = pickle.load(f)
            # A normal pnm file is a dict of lists (projects)
            if type(d) is dict:
                for name in d.keys():
                    # If dict item is a list, assume it's a valid project
                    if isinstance(d[name], list):
                        temp[name] = d[name]
                    else:
                        warnings.warn('File contents must be a dictionary, ' +
                                      'of lists, or a single list')
            else:
                if isinstance(d, list):  # If pickle contains a single list
                    temp[filename] = d
                else:
                    warnings.warn('File contents must be a dictionary, ' +
                                  'of lists, or a single list')

        # Now scan through temp dict to ensure valid types and names
        conflicts = set(temp.keys()).intersection(set(self.keys()))
        for name in list(temp.keys()):
            if name in conflicts:
                new_name = self._gen_name()
                warnings.warn('A project named ' + name + ' already exists, ' +
                              'renaming to ' + new_name)
                self[new_name] = temp[name]
            else:
                self[name] = temp[name]

    def _parse_filename(self, filename, ext='pnm'):
        p = Path(filename)
        p = p.resolve()
        # If extension not part of filename
        ext = ext.strip('.')
        if p.suffix != ('.' + ext):
            p = p.with_suffix(p.suffix + '.' + ext)
        return p

    def close_project(self, project):
        r"""
        Removes the specified Project from the Workspace

        This does not save the project, so any changes will be lost.
        """
        del self[project.name]

    def copy_project(self, project, new_name=None):
        r"""
        Make a copy of an existing Project

        Parameters
        ----------
        project : Project object
            The Project object to be copied

        name : string, optional
            A name for the new copy of the project.  If not supplied, then
            one will be generated (e.g. 'proj_02')

        Returns
        -------
        A handle to the new Project

        """
        new_sim = copy.deepcopy(project)
        new_sim._name = hex(id(new_sim))  # Assign temporary name
        new_sim.name = new_name
        return new_sim

    def new_project(self, name=None):
        r"""
        Creates a new empty Project object

        Parameters
        ----------
        name : string (optional)
            The unique name to give to the project.  If none is given, one
            will be automatically generated (e.g. 'sim_01`)

        Returns
        -------
        An empty project object, suitable for passing into a Network
        generator

        """
        sim = openpnm.utils.Project(name=name)
        return sim

    def _gen_name(self):
        r"""
        Generates a valid name for projects
        """
        n = [0]
        for item in self.keys():
            if item.startswith('sim_'):
                n.append(int(item.split('sim_')[1]))
        name = 'sim_'+str(max(n)+1).zfill(2)
        return name

    def _gen_ids(self, size):
        r"""
        Generates a sequence of integers of the given ``size``, starting at 1
        greater than the last produced value.

        The Workspace object keeps track of the most recent value, which
        persists until the current python session is restarted, so the
        returned array contains unique values for the given session.

        Parameters
        ----------
        size : int
            The number of values to generate.

        Returns
        -------
        A Numpy array of the specified size, containing integer values starting
        from the last used values.

        Notes
        -----
        When a new Workspace is created the
        """
        if not hasattr(self, '_next_id'):
            # If _next_id has not been set, then assign it
            self._next_id = 0
            # But check ids in any objects present first
            for proj in self.values():
                if len(proj) > 0:
                    if 'pore._id' in proj.network.keys():
                        Pmax = proj.network['pore._id'].max() + 1
                        Tmax = proj.network['throat._id'].max() + 1
                        self._next_id = max([Pmax, Tmax, self._next_id])
        ids = np.arange(self._next_id, self._next_id + size, dtype=np.int64)
        self._next_id += size
        return ids

    def __str__(self):
        s = []
        hr = '―'*78
        s.append(hr)
        s.append('OpenPNM Version ' + openpnm.__version__ + ' Workspace')
        s.append(hr)
        for item in self.values():
            s.append(' ' + item.name)
            s.append(item.__str__())
        return '\n'.join(s)
