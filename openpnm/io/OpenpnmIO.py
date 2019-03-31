import pickle
import warnings
import time
from flatdict import FlatDict
from openpnm.utils import NestedDict, sanitize_dict, Workspace
from openpnm.utils import logging
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)
ws = Workspace()


class OpenpnmIO(GenericIO):
    r"""
    This class possesses methods used for saving and loading OpenPNM
    workspaces, projects, and objects.

    """

    @classmethod
    def save_project(cls, project, filename=''):
        if filename == '':
            filename = project.name
        filename = cls._parse_filename(filename=filename, ext='pnm')

        # Save dictionary as pickle
        d = {project.name: project}
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    @classmethod
    def save_workspace(cls, filename=''):
        if filename == '':
            filename = 'workspace' + '_' + time.strftime('%Y%b%d_%H%M%p')
        filename = cls._parse_filename(filename=filename, ext='pnm')
        d = {}
        for sim in ws.values():
            d[sim.name] = sim
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    @classmethod
    def load_workspace(cls, filename, overwrite=False):

        fname = cls._parse_filename(filename=filename, ext='pnm')

        temp = {}  # Read file into temporary dict
        with open(fname, 'rb') as f:
            d = pickle.load(f)
            # A normal pnm file is a dict of lists (projects)
            if isinstance(d, dict):
                for name in d.keys():
                    # If dict item is a list, assume it's a valid project
                    if isinstance(d[name], list):
                        temp[name] = d[name]
                    else:
                        warnings.warn('File contents must be a dictionary ' +
                                      'of lists, or a single list')
        if overwrite:
            ws.clear()
        # Now scan through temp dict to ensure valid types and names
        conflicts = set(temp.keys()).intersection(set(ws.keys()))
        for name in list(temp.keys()):
            if name in conflicts:
                new_name = ws._gen_name()
                warnings.warn('A project named ' + name + ' already exists, ' +
                              'renaming to ' + new_name)
                ws[new_name] = temp[name]
            else:
                ws[name] = temp[name]
        return ws

    @classmethod
    def load_project(cls, filename):
        filename = cls._parse_filename(filename=filename, ext='pnm')
        temp = {}  # Read file into temporary dict
        with open(filename, 'rb') as f:
            d = pickle.load(f)
            # A normal pnm file is a dict of lists (projects)
            if type(d) is dict:
                cls.load_workspace(filename=filename, overwrite=False)
            elif isinstance(d, list):  # If pickle contains a single list
                temp[filename] = d
            else:
                warnings.warn('File contents must be a dictionary  ' +
                              'of lists, or a single list')

        # Now scan through temp dict to ensure valid types and names
        conflicts = set(temp.keys()).intersection(set(ws.keys()))
        for name in list(temp.keys()):
            if name in conflicts:
                new_name = ws._gen_name()
                warnings.warn('A project named ' + name + ' already exists, ' +
                              'renaming to ' + new_name)
                ws[new_name] = temp[name]
            else:
                ws[name] = temp[name]
        return ws
