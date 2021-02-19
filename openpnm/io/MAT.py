import scipy.io as spio
from flatdict import FlatDict
from openpnm.io import GenericIO, Dict
from openpnm.utils import sanitize_dict, logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


class MAT(GenericIO):
    r"""
    MAT files are a format used by Matlab

    Notes
    -----
    The 'mat' file must contain data formatted as follows:

    1. The file can contain either or both pore and throat data.

    2. The property names should be in the format of ``pore_volume`` or
    ``throat_surface_area``. In OpenPNM the first '_' will be replaced by
    a '.' to give ``'pore.volume'`` or ``'throat.surface_area'``.

    3. Boolean data represented as 1's and 0's will be converted to the
    Python boolean True and False.  These will become \'labels\' in
    OpenPNM.

    """

    @classmethod
    def save(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Using ``export_data`` instead.
        """
        cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, network, phases=[], filename=''):
        r"""
        Write Network to a Mat file for exporting to Matlab.

        Parameters
        ----------
        network : OpenPNM Network Object

        filename : string
            Desired file name, defaults to network name if not given

        phases : list of phase objects ([])
            Phases that have properties we want to write to file

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        network = network[0]
        # Write to file
        if filename == '':
            filename = project.name
        filename = cls._parse_filename(filename=filename, ext='mat')

        d = Dict.to_dict(network=network, phases=phases, interleave=True)
        d = FlatDict(d, delimiter='|')
        d = sanitize_dict(d)
        new_d = {}
        for key in list(d.keys()):
            new_key = key.replace('|', '_').replace('.', '_')
            new_d[new_key] = d.pop(key)

        spio.savemat(file_name=filename, mdict=new_d)

    @classmethod
    def load(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Using ``import_data`` instead.
        """
        return cls.import_data(*args, *kwargs)

    @classmethod
    def import_data(cls, filename, project=None):
        r"""
        Loads data from an appropriately formatted 'mat' file
        (i.e. MatLAB output).

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project object is supplied then one will be created and
            returned.

        Returns
        -------
        If no project object is supplied then one will be created and returned.

        """
        filename = cls._parse_filename(filename=filename, ext='mat')
        data = spio.loadmat(filename)
        # Reinsert the '.' separator into the array names
        for item in list(data.keys()):
            if item in ['__header__', '__version__', '__globals__']:
                data.pop(item)
                continue
            elif '_pore_' in item:
                path, prop = item.split('_pore_')
                new_key = path + '|pore.' + prop
            elif '_throat_' in item:
                path, prop = item.split('_throat_')
                new_key = path + '|throat.' + prop
            data[new_key] = data.pop(item)

        if project is None:
            project = ws.new_project()
        project = Dict.from_dict(data, project=project, delim='|')

        project = cls._convert_data(project)

        return project
