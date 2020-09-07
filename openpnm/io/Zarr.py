from zarr import DirectoryStore, ZipStore, group
from openpnm.io import GenericIO
from openpnm.utils import Workspace, Project, logging
ws = Workspace()


class Zarr(GenericIO):

    @classmethod
    def export_data(cls, objs, filename=None, storage='dir'):
        r"""
        export all data to a Zarr format

        Parameters
        ----------
        objs : list of OpenPNM objects
            The objects whose data should be added to the Zarr file
        filename : string
            The name of file
        storage : string
            The type of Zarr file to create.  Options are 'dir' for
            ``DirectoryStore`` (default) or 'zip' for ``ZipStore``.  For
            other options use the Zarr package directly.

        """
        if not isinstance(objs, list):
            objs = [objs]
        project = objs[0].project
        if filename is None:
            filename = project.name

        close = False
        if storage == 'dir':
            filename = cls._parse_filename(filename=filename, ext='zarr.dir')
            storage = DirectoryStore
        elif storage == 'zip':
            filename = cls._parse_filename(filename=filename, ext='zarr.zip')
            storage = ZipStore
            close = True
        else:
            raise Exception('To use more advanced formats use Zarr directly')

        store = storage(filename)
        root = group(store=store, overwrite=True)
        for obj in objs:
            item = root.create_group(obj.name, overwrite=True)
            item.update(obj)

        if close:
            store.close()
