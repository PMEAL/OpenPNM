import pickle
import openpnm
import time
import copy
import warnings
from pathlib import Path
from openpnm.core import logging
from openpnm.utils import SettingsDict
logger = logging.getLogger()


def singleton(cls):
    r'''
    This is based on PEP318
    '''
    instances = {}

    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance


@singleton
class Workspace(dict):

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
        if not isinstance(project, openpnm.core.Project):
            project = openpnm.core.Project(project, name=name)
        super().__setitem__(name, project)

    def _setloglevel(self, level):
        logger.setLevel(level)

    def _getloglevel(self):
        return 'Log level is currently set to: ' + str(logger.level)

    loglevel = property(fget=_getloglevel, fset=_setloglevel)

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
        Saves all the current *Projects* to a single 'pnm' file.

        Parameters
        ----------
        filename : string, optional
            If no filename is given, the a name is genrated using the current
            time and date. See Notes for more information on valid file names.

        See Also
        --------
        save_project

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
            filename = 'workspace' + '_' + time.strftime('%Y%b%d_%H%M%p')
        filename = self._parse_filename(filename=filename, ext='pnm')
        d = {}
        for sim in self.values():
            d[sim.name] = sim
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    def load_workspace(self, filename):
        r"""
        Loads a saved OpenPNM session from a 'pnm' file.  If the 'pnm' file
        contains multiple *Projects*, they will all be loaded.  Any *Projects*
        present in the current *Workspace* will be deleted.

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
        Save given project to a 'pnm' file, including all of associated
        objects including algorithms.

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
        Loads a *Project* from the specified 'pnm' file and adds it to the
        *Workspace*.  This will *not* delete any existing *Projects* in the
        *Workspace* and will rename any *Projects* being loaded if necessary.

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
        """
        del self[project.name]

    def copy_project(self, project, new_name=None):
        r"""
        Make a copy of an existing project

        Parameters
        ----------
        project : Project object
            The Project object to be copied

        name : string, optional
            A name for the new copy of the project.  If not supplied, then
            one will be generated (e.g. 'proj_02')

        Returns
        -------
        The new Project object

        """
        new_sim = copy.deepcopy(project)
        new_sim._name = hex(id(new_sim))  # Assign temporary name
        new_sim.name = new_name
        return new_sim

    def new_project(self, name=None):
        r"""
        Creates a new, empty project object

        Parameters
        ----------
        name : string (optional)
            The unique name to give to the project.  If none is given, one
            will be automatically generated (e.g. 'sim_01`)

        Returns
        -------
        An empty project object, sui+table for passing into a Network
        generator

        """
        sim = openpnm.core.Project(name=name)
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
