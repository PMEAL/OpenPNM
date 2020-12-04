import openpnm
import numpy as np
from openpnm.utils import SettingsDict, logging
logger = logging.getLogger(__name__)


class SettingsDict(SettingsDict):
    def __setitem__(self, key, value):
        if key == 'loglevel':
            logger = logging.getLogger()
            logger.setLevel(value)
        super().__setitem__(key, value)


class Workspace(dict):
    r"""
    The Workspace object provides the highest level of adminstrative
    control over active OpenPNM sessions.

    It is a ``dictionary`` that stores a list of all open Projects by name.

    This class is a
    `singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ so that
    whenever and wherever a Workspace is instantiated, the same instance is
    obtained.  This allows it to maintain a definitive record of all open
    Projects.

    See Also
    --------
    Project

    Notes
    -----
    The Workspace object contains a variety of functions that one might expect
    from the 'file-menu' in a typical GUI.

    """

    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Workspace.__instance__ is None:
            Workspace.__instance__ = dict.__new__(cls)
        return Workspace.__instance__

    def __init__(self):
        super().__init__()
        self._projects = {}
        self.settings = SettingsDict()
        self.settings['loglevel'] = 30

    def __setitem__(self, name, project):
        if name is None:
            name = self._gen_name()
        if name in self.keys():
            # To ensure renaming a Project to its current name is allowed
            if self.__instance__[name] is not project:
                raise Exception(f"A project named {name} already exists")
        if project in self.values():
            self.pop(project.name, None)
        if not isinstance(project, openpnm.utils.Project):
            project = openpnm.utils.Project(project, name=name)
        super().__setitem__(name, project)

    def copy(self):
        r"""
        """
        raise Exception('Cannot copy Workspace, only one can exist at a time')

    def _create_console_handles(self, project):
        r"""
        Adds all objects in the given project to the console as variables
        with handle names taken from each object's name.
        """
        import __main__
        for item in project:
            __main__.__dict__[item.name] = item

    @property
    def version(self):
        return openpnm.__version__

    def save_workspace(self, filename=None):
        r"""
        Save all projects in the current workspace as a single file

        Parameters
        ----------
        filename : str
            The filename to use when saving.  If not provided, the present
            date and time are used.

        Notes
        -----
        The file is actually zip archive containing ``pnm`` files, one for
        each project in the workspace. This archive can be extracted and each
        ``pnm`` file can be loaded manually using ``load_project`` or the
        ``openpnm.io.PNM`` class.
        """
        from datetime import datetime
        if filename is None:
            dt = datetime.now()
            filename = dt.strftime("%Y_%m_%d_%H_%M_%S")
        from zipfile import ZipFile
        with ZipFile(filename + '.wrk', 'w') as z:
            for prj in self.values():
                prj.save_project()
                z.write(prj.name + '.pnm')

    def load_workspace(self, filename):
        r"""
        Load project(s) from a saved workspace into current workspace

        Parameters
        ----------
        filename : str or path object
            The filename containing the saved workspace

        Notes
        -----
        ??

        """
        from zipfile import ZipFile
        with ZipFile(filename, 'r') as z:
            logger.info('Loading projects contained in ' + filename)
            files = z.filelist
            for f in files:
                self.load_project(f.orig_filename)

    def save_project(self, project, filename=None):
        r"""
        Saves given Project to a ``pnm`` file

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
        from openpnm.io import PNM
        PNM.save_project(project=project, filename=filename)

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
        from openpnm.io import PNM
        proj = PNM.load_project(filename=filename)
        return proj

    def close_project(self, project):
        r"""
        Removes the specified Project from the Workspace

        This does not save the project, so any changes will be lost.
        """
        del self[project.name]

    def copy_project(self, project, name=None):
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
        proj : list
            A handle to the new Project

        """
        proj = project.copy(name)
        return proj

    def new_project(self, name=None):
        r"""
        Creates a new empty Project object

        Parameters
        ----------
        name : string (optional)
            The unique name to give to the project.  If none is given, one
            will be automatically generated (e.g. 'proj_01`)

        Returns
        -------
        proj : list
            An empty Project object, suitable for passing into a Network
            generator

        """
        proj = openpnm.utils.Project(name=name)
        return proj

    def _gen_name(self):
        r"""
        Generates a valid name for projects
        """
        n = [0]
        for item in self.keys():
            if item.startswith('proj_'):
                n.append(int(item.split('proj_')[1]))
        name = 'proj_'+str(max(n)+1).zfill(2)
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
