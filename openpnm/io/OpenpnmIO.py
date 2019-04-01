import pickle
import time
from openpnm.utils import Workspace
from openpnm.utils import logging
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)
ws = Workspace()


class OpenpnmIO(GenericIO):
    r"""
    This class possesses methods used for saving and loading OpenPNM
    workspaces, projects, and objects.

    Notes
    -----
    The methods in this class use the ``pickle`` module from the standard
    library, which have known security issues.  Do not open '.pnm' files
    from untrusted sources.

    """

    @classmethod
    def save_objects(cls, objs):
        if not isinstance(objs, list):
            objs = [objs]
        for item in objs:
            fname = cls._parse_filename(filename=item.name,
                                        ext=item.settings['prefix'])
            with open(fname, 'wb') as f:
                pickle.dump({item.name: item}, f)

    @classmethod
    def load_object(cls, filename, project):
        p = cls._parse_filename(filename)
        with open(p, 'rb') as f:
            d = pickle.load(f)
        obj = project._new_object(objtype=p.suffix.strip('.'),
                                  name=p.name.split('.')[0])
        obj.update(d)

    @classmethod
    def load(cls, filename):
        fname = cls._parse_filename(filename=filename, ext='pnm')
        with open(fname, 'rb') as f:
            d = pickle.load(f)
        if isinstance(d, dict):
            ws = cls.load_workspace(filename=filename, overwrite=True)
            return ws
        elif isinstance(d, list):
            proj = cls.load_project(filename=filename)
            return proj

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
                        raise Exception('File does not contain a valid ' +
                                        'OpenPNM Workspace')

        if overwrite:
            ws.clear()
        # Now scan through temp dict to ensure valid types and names
        conflicts = set(temp.keys()).intersection(set(ws.keys()))
        for name in list(temp.keys()):
            if name in conflicts:
                new_name = ws._gen_name()
                logger.warning('A project named ' + name + ' already exists,' +
                               ' renaming to ' + new_name)
                ws[new_name] = temp[name]
            else:
                ws[name] = temp[name]
        return ws

    @classmethod
    def load_project(cls, filename):
        filename = cls._parse_filename(filename=filename, ext='pnm')
        projname = filename.name.split('.')[0]
        with open(filename, 'rb') as f:
            d = pickle.load(f)
            # A normal pnm file is a dict of lists (projects), so check first
            if type(d) is dict:
                if len(d) == 1:  # If only one project in dict
                    # Store existing keys found in Workspace
                    projects = set(ws.keys())
                    # Call load_workspace to handle file
                    temp = cls.load_workspace(filename=filename,
                                              overwrite=False)
                    # Find new addition to Workspace and return it
                    proj = list(set(temp.keys()).difference(projects))[0]
                    # Return Project handle to user and exit
                    return ws[proj]
                else:
                    raise Exception(filename.name + ' contains multiple ' +
                                    'projects, use load_workspace instead')
            # If pickle contains a single list
            elif isinstance(d, list):
                if projname not in ws.keys():
                    ws[projname] = d
                    return ws[projname]
                else:
                    newname = ws._gen_name()
                    logger.warning('Project named ' + projname +
                                   ' already present in Workspace,' +
                                   ' renaming to ' + newname)
                    ws[newname] = d
                    return ws[newname]
            else:
                raise Exception('File contents are not understood')
