import zarr
import numpy as np
import importlib
from datetime import datetime
from openpnm.utils import Workspace, Project
from openpnm.utils import logging
from openpnm.io import GenericIO
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

        with zarr.ZipStore(filename, mode='w') as store:
            root = zarr.group(store=store)

            # root.attrs['version'] = op.__version__
            root.attrs['date saved'] = datetime.today().strftime("%Y %h %d %H:%M:%S")
            root.attrs['comments'] = project.comments
            for obj in project:
                item = root.create_group(obj.name)
                # Store data
                item.update(obj)
                # Store settings dict as metadata
                item.attrs['settings'] = obj.settings
                # Store models dict as metadata
                if hasattr(obj, 'models'):
                    temp = obj.models.copy()
                    for m in temp.keys():
                        a = temp[m].pop('model', None)
                        if a is not None:
                            temp[m]['model'] = a.__module__ + '|' + \
                                a.__code__.co_name
                    item.attrs['models'] = temp
                item.attrs['class'] = str(obj.__class__)

    @classmethod
    def load_project(cls, filename):
        loglevel = ws.settings['loglevel']
        ws.settings['loglevel'] = 50
        with zarr.ZipStore(filename, mode='r') as store:
            root = zarr.group(store=store)
            proj = Project()
            for name in root.keys():
                if 'network' in root[name].attrs['class']:
                    proj, obj = create_obj(root, name, proj)
            for name in root.keys():
                if 'network' not in root[name].attrs['class']:
                    proj, obj = create_obj(root, name, proj)

        ws.settings['loglevel'] = loglevel
        return proj


def create_obj(root, name, proj):
    import openpnm as op
    # regenerate object as same class
    mro = root[name].attrs['class']
    mro = mro.split("'")[1]
    mro = mro.split('.')
    mod = getattr(op, mro[1])
    c = [i for i in mod.__dir__() if i.startswith('Generic')][0]
    c = mro[-1]
    clss = getattr(mod, c)
    obj = clss(project=proj)
    obj._name = name
    # Add data to obj
    for item in root[name]:
        obj.update({item: np.array(root[name][item])})
    # Add settings to obj
    obj.settings.update(root[name].attrs['settings'])
    # Add models to obj
    if hasattr(obj, 'models'):
        obj.models.update(root[name].attrs['models'])
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
    return proj, obj
