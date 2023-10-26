import logging
import pickle
import re
import sys
from datetime import datetime
from uuid import uuid4

from openpnm.utils import SettingsAttr

logger = logging.getLogger("openpnm")

__all__ = ['Workspace']


class WorkspaceSettings(SettingsAttr):
    r"""

    Parameters
    ----------
    default_solver : str
        The solver to use by default, if user does not specify one explicitly.
        The default values is PardisoSpsolve, but a good option is ScipySpsolve
        if the Pardiso is causing problems.
    loglevel : int
        Sets the threshold for the severity of logger message which appear.
        Ranges are as follows:

        ======= ==============================================================
        Level   Description
        ======= ==============================================================
        10      DEBUG: Detailed information, typically of interest only when
                diagnosing problems.
        20      INFO: Confirmation that things are working as expected.
        30      WARNING: An indication that something unexpected happened, or
                indicative of some problem in the near future (e.g. ‘disk
                space low’). The software is still working as expected.
        40      ERROR: Due to a more serious problem, the software has not
                been able to perform some function.
        50      CRITICAL: A serious error, indicating that the program itself
                may be unable to continue running.
        ======= ==============================================================
    """
    # Pardiso requires MKL, which is not available on new Apple chips
    if sys.platform == 'darwin':
        default_solver = 'ScipySpsolve'
    else:
        try:
            import pypardiso
            default_solver = 'PardisoSpsolve'
        except ImportError:
            default_solver = 'ScipySpsolve'
            msg = (
                'PARDISO solver not installed, run `pip install pypardiso`. '
                'Otherwise, simulations will be slow. Apple M chips not supported.'
            )
            logger.error(msg)

    @property
    def loglevel(self):
        return self._loglevel

    @loglevel.setter
    def loglevel(self, value):
        if isinstance(value, str):
            options = {
                "TRACE": 5,
                "DEBUG": 10,
                "INFO": 20,
                "SUCESS": 25,
                "WARNING": 30,
                "ERROR": 40,
                "CRITICAL": 50
            }
            value = options[value]
        self._loglevel = value
        logger.setLevel(value)


class Workspace(dict):
    r"""
    The Workspace object provides the highest level of adminstrative
    control over active OpenPNM sessions.

    It is a ``dictionary`` that stores a list of all open Projects by name.

    This class is a
    `singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ so that
    whenever and wherever a Workspace is instantiated, the same instance is
    obtained. This allows it to maintain a definitive record of all open
    Projects.

    See Also
    --------
    Project

    Notes
    -----
    The Workspace object contains a variety of functions that one might
    expect from the 'file-menu' in a typical GUI.

    """

    __instance__ = None

    # This __new__ method makes the Workspace a Singleton
    def __new__(cls, *args, **kwargs):
        if Workspace.__instance__ is None:
            Workspace.__instance__ = dict.__new__(cls)
        return Workspace.__instance__

    def __init__(self):
        super().__init__()
        self.settings = WorkspaceSettings()
        self.settings.loglevel = 30

    def copy(self):
        """Brief explanation of 'copy'"""
        raise Exception('Cannot copy Workspace, only one can exist at a time')

    def _create_console_handles(self, project):
        r"""
        Adds all objects in the given project to the console as variables
        with handle names taken from each object's name.
        """
        import __main__
        for item in project:
            __main__.__dict__[item.name] = item

    def _validate_name(self, name=None):
        r"""
        Generates a valid name for projects
        """
        if name in [None, '']:
            name = 'proj_01'  # Give basic name, then let rest of func fix it
        if name in self.keys():  # If proposed name is taken, increment it
            if not re.search(r'_\d+$', name):  # If name does not end with _##
                name = name + '_01'
            prefix, count = name.rsplit('_', 1)
            n = [0]
            for item in list(self.keys()):
                if item.startswith(prefix+'_'):
                    n.append(int(item.split(prefix+'_')[1]))
            name = prefix+'_'+str(max(n)+1).zfill(2)
        return name

    def save_workspace(self, filename=None):
        r"""
        Saves all projects in the current workspace as a single file

        Parameters
        ----------
        filename : str
            The filename to use when saving. If not provided, the present
            date and time are used.

        Notes
        -----
        The file is actually zip archive containing ``pnm`` files, one for
        each project in the workspace. This archive can be extracted and each
        ``pnm`` file can be loaded manually using ``load_project`` or the
        ``openpnm.io.PNM`` class.

        """
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
        Loads project(s) from a saved workspace into current workspace

        Parameters
        ----------
        filename : str or Path
            The filename containing the saved workspace

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
        project : Project
            The project to save.
        filename : str, optional
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
        if filename is None:
            dt = datetime.now()
            filename = dt.strftime("%Y_%m_%d_%H_%M_%S")
        with open(filename.split('.')[-1]+'.pnm', 'wb') as f:
            pickle.dump(project, f)

    def load_project(self, filename):
        r"""
        Loads a Project from the specified 'pnm' file

        The loaded project is added to the Workspace. This will *not* delete
        any existing Projects in the Workspace and will rename any Projects
        being loaded if necessary.

        Parameters
        ----------
        filename : str or Path
            The name of the file to open. See Notes for more information.

        See Also
        --------
        load_workspace

        Notes
        -----
        The filename can be a string such as 'saved_file.pnm'. The string can
        include absolute path such as 'C:\networks\saved_file.pnm', or can
        be a relative path such as '..\..\saved_file.pnm', which will look
        2 directories above the current working directory.  Can also be a
        path object object such as that produced by ``pathlib`` or
        ``os.path`` in the Python standard library.

        """
        with open(filename, 'rb') as f:
            proj = pickle.load(f)
            proj.settings.uuid = str(uuid4())
            proj.name = self._validate_name(proj.name)
            self[proj.name] = proj
        return proj

    def close_project(self, project):
        r"""
        Removes the specified Project from the Workspace

        This does not save the project, so any changes will be lost.

        Parameters
        ----------
        project : Project
            The Project object to be copied

        """
        _ = self.pop(project.name, None)

    def copy_project(self, project, name=None):
        r"""
        Makes a copy of an existing Project

        Parameters
        ----------
        project : Project
            The Project object to be copied
        name : str, optional
            A name for the new copy of the project. If not supplied, then
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
        name : str, optional
            The unique name to give to the project. If none is given, one
            will be automatically generated (e.g. 'proj_01`)

        Returns
        -------
        proj : list
            An empty Project object, suitable for passing into a Network
            generator

        """
        from openpnm.utils import Project
        proj = Project(name=name)
        return proj

    def __str__(self):  # pragma: no cover
        hr = '―'*78
        s = hr + '\n'
        for item in self.values():
            s += item.name + '\n'
            s += item.__str__() + '\n'
        return s
