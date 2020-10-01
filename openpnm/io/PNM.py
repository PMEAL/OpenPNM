import json
import json_tricks as jsont
import numpy as np
import importlib
from datetime import datetime
from openpnm.utils import Workspace, Project
from openpnm.utils import logging
from openpnm.io import GenericIO
from h5py import File as hdfFile
logger = logging.getLogger(__name__)
ws = Workspace()


class PNM(GenericIO):
    r"""
    This is the official way to save and load OpenPNM projects


    """

    @classmethod
    def save_project(cls, project, filename=None):
        if filename is None:
            filename = project.name + '.pnm'

        # Make a directory using the given file name
        f = cls._parse_filename(filename, 'pnm')
        with hdfFile(f, mode='w') as root:
            # root = hdfFile(f, mode='w')
            root.attrs['version'] = ws.version
            date = datetime.today().strftime("%Y %h %d %H:%M:%S")
            root.attrs['date saved'] = date
            root.attrs['name'] = project.name
            # root.attrs['comments'] = project.comments
            for obj in project:
                item = root.create_group(obj.name)
                for arr in obj.keys():  # Store data
                    try:
                        item.create_dataset(name=arr, data=obj[arr],
                                            shape=obj[arr].shape,
                                            compression="gzip")
                    except TypeError:  # Deal with 'object' arrays
                        logger.warning(arr + ' is being converted to a string')
                        b = jsont.dumps(obj[arr])
                        c = b.encode()
                        d = np.void(c)
                        item.create_dataset(name=arr, data=d)
                # Store settings dict as metadata
                item.attrs['settings'] = json.dumps(obj.settings)
                # Store models dict as metadata
                if hasattr(obj, 'models'):
                    obj_models = {}
                    for model in obj.models.keys():
                        temp = {k: v for k, v in obj.models[model].items()
                                if k != 'model'}
                        if 'model' in obj.models[model].keys():
                            a = obj.models[model]['model']
                            temp['model'] = a.__module__ + '|' + \
                                a.__code__.co_name
                        obj_models[model] = temp
                    item.attrs['models'] = json.dumps(obj_models)
                item.attrs['class'] = str(obj.__class__)

    @classmethod
    def load_project(cls, filename):
        loglevel = ws.settings['loglevel']
        ws.settings['loglevel'] = 50
        f = cls._parse_filename(filename, 'pnm')
        with hdfFile(f, mode='r') as root:
            try:  # Create an empty project with old name
                proj = Project(name=root.attrs['name'])
                print('Loading ' + proj.name)
            except Exception:  # Generate a new name if collision occurs
                proj = Project()
                print('A project named ' + root.attrs['name']
                      + ' already exists, renaming to ' + proj.name)
            print('Created using OpenPNM version ' + root.attrs['version'])
            print('Saved on ' + root.attrs['date saved'])
            for name in root.keys():
                if 'network' in root[name].attrs['class']:
                    proj, obj = create_obj(root, name, proj)
            for name in root.keys():
                if 'network' not in root[name].attrs['class']:
                    proj, obj = create_obj(root, name, proj)
            ws.settings['loglevel'] = loglevel
        return proj


def create_obj(root, name, proj):
    r"""
    Reproduces an OpenPNM object, given the hdf5 file and name
    """
    import openpnm as op
    # regenerate object as same class
    mro = root[name].attrs['class']
    mro = mro.split("'")[1]
    mro = mro.split('.')
    mod = importlib.import_module('.'.join(mro[:-1]))
    clss = getattr(mod, mro[-1])
    obj = clss(project=proj, settings={'freeze_models': True})
    obj._name = name
    # Add data to obj
    for arr in root[name].keys():
        a = np.array(root[name][arr])
        if str(a.dtype).startswith("|V"):
            logger.warning(arr + ' is being converted from string')
            b = np.string_(a)
            c = b.astype(str)
            a = jsont.loads(c)
        obj.update({arr: a})
    # Add settings to obj
    obj.settings.update(json.loads(root[name].attrs['settings']))
    # Add models to obj
    if hasattr(obj, 'models'):
        obj.models.update(json.loads(root[name].attrs['models']))
        for m in obj.models.keys():
            md, fn = obj.models[m]['model'].split('|')
            try:
                md = importlib.import_module(md)
                try:
                    obj.models[m]['model'] = getattr(md, fn)
                except AttributeError:
                    print('Warning: the function \"' + fn
                          + '\" could not be loaded, adding \"blank\" instead')
                    obj.models[m]['model'] = op.models.misc.basic_math.blank
            except ModuleNotFoundError:
                print('Warning: the module \"' + md
                      + '\" could not be found, adding \"blank\" instead')
                obj.models[m]['model'] = op.models.misc.basic_math.blank
    del obj.settings['freeze_models']
    return proj, obj
