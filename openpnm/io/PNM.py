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
                found_attrs = set(obj.__dict__.keys())
                known_attrs = set(['settings', '_models_dict',
                                   '_am', '_im',
                                   '_spacing', '_shape'])
                foreign_attrs = found_attrs.difference(known_attrs)
                if len(foreign_attrs) > 0:
                    line_break = f"\n{'':13}"
                    logger.critical(f"{obj.name} has the following attributes that will"
                                    + f" not be saved: {[i for i in foreign_attrs]}"
                                    + f"{line_break}Consider using Pickle instead")
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
                    model = None
                    for model in obj.models.keys():
                        temp = {k: v for k, v in obj.models[model].items()
                                if k != 'model'}
                        if 'model' in obj.models[model].keys():
                            a = obj.models[model]['model']
                            temp['model'] = a.__module__ + '|' + \
                                a.__code__.co_name
                        obj_models[model] = temp
                    try:
                        item.attrs['models'] = json.dumps(obj_models)
                    except TypeError:
                        logger.critical(f'The model {model} and its parameters'
                                        + ' could not be written to file.')
                item.attrs['class'] = str(obj.__class__)

    @classmethod
    def load_project(cls, filename):
        f = cls._parse_filename(filename, 'pnm')
        with hdfFile(f, mode='r') as root:
            logger.info('Loading project from file ' + f.name)
            try:  # Create an empty project with old name
                proj = Project(name=root.attrs['name'])
                logger.info('Loading ' + proj.name)
            except Exception:  # Generate a new name if collision occurs
                proj = Project()
                logger.warning('A project named ' + root.attrs['name']
                               + ' already exists, renaming to ' + proj.name)
            logger.info('Created using OpenPNM version ' + root.attrs['version'])
            logger.info('Saved on ' + root.attrs['date saved'])
            loglevel = ws.settings['loglevel']
            ws.settings['loglevel'] = 50
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
    obj = clss.__new__(cls=clss)
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
                    logger.warning(f"The function {fn} could not be loaded, adding"
                                   + " 'blank' instead")
                    obj.models[m]['model'] = op.models.misc.basic_math.blank
            except ModuleNotFoundError:
                logger.warning(f"The module {md} could not be loaded, adding"
                               + " 'blank' instead")
                obj.models[m]['model'] = op.models.misc.basic_math.blank
    proj.append(obj)
    return proj, obj
