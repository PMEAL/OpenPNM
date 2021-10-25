import time
import pickle
from openpnm.io import GenericIO
from openpnm.utils import logging
from openpnm.utils import Workspace, Project
logger = logging.getLogger(__name__)
ws = Workspace()


class Pickle(GenericIO):
    r"""
    This class contains methods used for saving and loading OpenPNM Workspaces,
    Projects, and objects as Pickles.

    Notes
    -----
    The methods in this class use the ``pickle`` module from the standard
    library.  Aside from the security issues, these files can only be loaded
    by the exact same OpenPNM version used to save them.  They are meant for
    **temporary storage**.

    """

    @classmethod
    def save_object_to_file(cls, objs):
        r"""
        Saves an OpenPNM object or list of objects to a file of set of files

        Parameters
        ----------
        objs : OpenPNM Base object or list of objects
            The object(s) to be saved

        """
        if not isinstance(objs, list):
            objs = [objs]
        for item in objs:
            fname = cls._parse_filename(filename=item.name,
                                        ext=item.settings['prefix'])
            with open(fname, 'wb') as f:
                pickle.dump({item.name: item}, f)

    @classmethod
    def load_object_from_file(cls, filename, project=None):
        r"""
        Loads an OpenPNM object from a file on disk

        Parameters
        ----------
        filename : string or path object
            The name of the file containing the object to open. Can be a
            relative or absolute path, and can be a string or path object such
            as that produced by ``pathlib``.
        project : OpenPNM Project object
            If not provided one will be created and returned by this function,
            otherwise the loaded object will be added to the given ``project``.

        Returns
        -------
        project : OpenPNM Project
            If no Project object is specified then one is created. A handle to
            the Project is returned.

        """
        if project is None:
            project = Project()
        p = cls._parse_filename(filename)
        with open(p, 'rb') as f:
            d = pickle.load(f)
        obj = project._new_object(objtype=p.suffix.strip('.'),
                                  name=p.name.split('.')[0])
        obj.update(d)
        return project

    @classmethod
    def save_project(cls, project, filename=''):
        r"""
        Save an OpenPNM Project to a file on disk

        Parameters
        ----------
        project : OpenPNM Project
            The project to save
        filename : string
            The filename to save the file
        """
        if filename == '':
            filename = project.name
        filename = cls._parse_filename(filename=filename, ext='pkl')

        # Save dictionary as pickle
        d = {project.name: project}
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    @classmethod
    def save_workspace(cls, filename=''):
        r"""
        Save the current Workspace to a file on disk

        Parameters
        ----------
        filename : string
            The filename to save the file
        """
        if filename == '':
            filename = 'workspace' + '_' + time.strftime('%Y%b%d_%H%M%p')
        filename = cls._parse_filename(filename=filename, ext='pkl')
        # Create a normal dict to store objects to prevent name errors upon
        # reopening
        d = {}
        for sim in ws.values():
            d[sim.name] = sim
        with open(filename, 'wb') as f:
            pickle.dump(d, f)

    @classmethod
    def load_workspace(cls, filename, overwrite=False):
        r"""
        Load a saved Workspace into the current one

        Parameters
        ----------
        filename : string or path object
            The name of the file to load
        overwrite : boolean
            A flag to indicate if the current Workspace should be
            overwritten when loading the new one.  The default is ``False``,
            meaning the loaded file will be added to the existing data.  Note
            that in this case Project names may clash, in which case the
            newly loaded Projects are given new names.

        Returns
        -------
        workspace : OpenPNM Workspace Object
            A handle to the Workspace, with the newly loaded Projects added
        """
        fname = cls._parse_filename(filename=filename, ext='pkl')
        temp = {}  # Read file into temporary dict
        if overwrite:
            ws.clear()
        with open(fname, 'rb') as f:
            d = pickle.load(f)
            # A normal pnm file is a dict of lists (projects)
            if isinstance(d, dict):
                for name in d.keys():
                    # If dict item is a list, assume it's a valid project
                    if isinstance(d[name], list):
                        temp[name] = d[name]
                    else:
                        raise Exception('File does not contain a valid '
                                        + 'OpenPNM Workspace')
        # Now scan through temp dict to ensure valid types and names
        conflicts = set(temp.keys()).intersection(set(ws.keys()))
        for name in list(temp.keys()):
            if name in conflicts:
                new_name = ws._gen_name()
                logger.warning('A project named ' + name + ' already exists,'
                               + ' renaming to ' + new_name)
                ws[new_name] = temp[name]
            else:
                ws[name] = temp[name]
        return ws

    @classmethod
    def load_project(cls, filename):
        r"""
        Load a saved Project file into the current Workspace

        Parameters
        ----------
        filename : string or path object
            The name of the file to load

        Returns
        -------
        project : OpenPNM Project
            A handle to the loaded Project is returned.
        """
        filename = cls._parse_filename(filename=filename, ext='pkl')
        projname = filename.name.split('.')[0]
        with open(filename, 'rb') as f:
            d = pickle.load(f)
            # A normal pnm file is a dict of lists (projects), so check first
            if isinstance(d, dict):
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
                raise Exception(filename.name + ' contains multiple'
                                + ' projects, use load_workspace instead')
            # If pickle contains a single list
            if isinstance(d, list):
                if projname not in ws.keys():
                    ws[projname] = d
                    return ws[projname]
                newname = ws._gen_name()
                logger.warning('Project named ' + projname
                               + ' already present in Workspace,'
                               + ' renaming to ' + newname)
                ws[newname] = d
                return ws[newname]
            # Otherwise, raise exception
            raise Exception('File contents are not understood')
