import zarr
import numpy as np
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
                    for m in temp:
                        temp[m].pop('model', None)
                    item.attrs['models'] = temp
                item.attrs['class'] = str(obj.__class__)

    @classmethod
    def load_project(cls, filename):
        loglevel = ws.settings['loglevel']
        ws.settings['loglevel'] =  50
        with zarr.ZipStore('example.pnm', mode='r') as store:
            root = zarr.group(store=store)
            proj = Project()
            for name in root.keys():
                if 'network' in root[name].attrs['class']:
                    proj, obj = create_obj(root, name, proj)
                    obj.settings.update(root[name].attrs['settings'])
                    obj.models.update(root[name].attrs['models'])
            for name in root.keys():
                if 'network' not in root[name].attrs['class']:
                    proj, obj = create_obj(root, name, proj)
                    obj.settings.update(root[name].attrs['settings'])
                    if hasattr(obj, 'models'):
                        obj.models.update(root[name].attrs['models'])
        ws.settings['loglevel'] =  loglevel
        return proj


def create_obj(root, name, proj):
    import openpnm as op
    mro = root[name].attrs['class']
    mro = mro.split("'")[1]
    mro = mro.split('.')
    mod = getattr(op, mro[1])
    c = [i for i in mod.__dir__() if i.startswith('Generic')][0]
    c = mro[-1]
    clss = getattr(mod, c)
    obj = clss(project=proj)
    obj._name = name
    for item in root[name]:
        obj.update({item: np.array(root[name][item])})
    return proj, obj
