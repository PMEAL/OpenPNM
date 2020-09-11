import zarr
import numpy as np
import importlib
from datetime import datetime
from openpnm.utils import Workspace, Project
from openpnm.utils import logging
from openpnm.io import GenericIO, XDMF
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

        store = zarr.DirectoryStore(filename)
        root = zarr.group(store=store, overwrite=True)
        root.attrs['version'] = ws.version
        date = datetime.today().strftime("%Y %h %d %H:%M:%S")
        root.attrs['date saved'] = date
        # root.attrs['comments'] = project.comments
        for obj in project:
            item = root.create_group(obj.name, overwrite=True)
            # Store data
            # item.update(obj)
            for arr in obj.keys():
                item.create_dataset(name=arr, data=obj[arr],
                                    shape=obj[arr].shape,
                                      compressor=None)
            # Store settings dict as metadata
            item.attrs['settings'] = obj.settings
            # Store models dict as metadata
            if hasattr(obj, 'models'):
                obj_models = {}
                for model in obj.models.keys():
                    temp = {k: v for k, v in obj.models[model].items() if k != 'model'}
                    if 'model' in obj.models[model].keys():
                        a = obj.models[model]['model']
                        temp['model'] = a.__module__ + '|' + \
                            a.__code__.co_name
                    obj_models[model] = temp
                item.attrs['models'] = obj_models
            item.attrs['class'] = str(obj.__class__)
        XDMF.save(network=project.network,
                  phases=list(project.phases().values()),
                  filename=filename + "/" + 'for_paraview.xmf')

    @classmethod
    def load_project(cls, filename):
        loglevel = ws.settings['loglevel']
        ws.settings['loglevel'] = 50
        store = zarr.DirectoryStore(filename)
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
    obj = clss(project=proj, settings={'freeze_models': True})
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
    del obj.settings['freeze_models']
    return proj, obj
