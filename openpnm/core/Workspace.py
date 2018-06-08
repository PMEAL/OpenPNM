import pickle
import openpnm
import time
import copy
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
        """
        if filename == '':
            filename = 'workspace' + '_' + time.strftime('%Y%b%d_%H%M%p')
        else:
            filename = filename.rsplit('.pnm', 1)[0]
        d = {}
        for sim in self.values():
            d[sim.name] = sim
        pickle.dump(d, open(filename + '.pnm', 'wb'))

    def load_workspace(self, filename):
        r"""
        """
        self.clear()
        self.load_project(filename=filename)

    def save_project(self, project, filename=''):
        r"""
        Save given project to a 'pnm' file, including all of associated
        objects, but not Algorithms.

        Parameters
        ----------
        project : OpenPNM Project
            The project to save

        filename : string, optional
            If no filename is given, the given project name is used
        """
        if filename == '':
            filename = project.name
        else:
            filename = filename.rsplit('.pnm', 1)[0]

        # Save dictionary as pickle
        d = {project.name: project}
        pickle.dump(d, open(filename + '.pnm', 'wb'))

    def load_project(self, filename):
        r"""
        Loads a project from the specified 'pnm' file and adds it
        to the Workspace

        Parameters
        ----------
        filename : string
            The name of the file containing the Network project to load
        """
        filename = filename.rsplit('.pnm', 1)[0]
        d = pickle.load(open(filename + '.pnm', 'rb'))
        if type(d) is dict:
            for name in d.keys():
                self[name] = d[name]
        elif type(d) is list:
            self[filename] = d
        else:
            raise Exception('Can only load a dict of lists or a list')

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
