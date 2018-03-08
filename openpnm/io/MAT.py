import scipy as sp
import scipy.io as spio
from openpnm.core import logging, Workspace
from openpnm.io import GenericIO, Dict
from openpnm.utils import FlatDict, sanitize_dict
logger = logging.getLogger(__name__)
ws = Workspace()


class MAT(GenericIO):
    r"""
    Class for reading and writing OpenPNM data to a Matlab 'mat' file

    Notes
    -----
    The 'mat' file must contain data formatted as follows:

    1. The file can contain either or both pore and throat data.

    2. The property names should be in the format of ``pore_volume`` or
    ``throat_surface_area`. In OpenPNM the first \'_\' will be replaced by
    a \'.\' to give \'pore.volume\' or \'throat.surface_area\'.

    3. Boolean data represented as 1's and 0's will be converted to the
    Python boolean True and False.  These will become \'labels\' in
    OpenPNM.
    """

    @classmethod
    def save(cls, network, phases=[], filename='', delim='_'):
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
        d = FlatDict(d, delimiter='.')
        d = sanitize_dict(d)

        spio.savemat(file_name=filename, mdict=d)

    @classmethod
    def load(cls, filename, project=None):
        r"""
        Loads data onto the given network from an appropriately formatted
        'mat' file (i.e. MatLAB output).

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
        If no Network object is supplied then one will be created and returned.

        If return_geometry is True, then a tuple is returned containing both
        the network and a geometry object.

        """
        filename = cls._parse_filename(filename=filename, ext='mat')
        data = spio.loadmat(filename)
        # Deal with pore coords and throat conns specially
        for item in list(data.keys()):
            if item in ['__header__', '__version__', '__globals__']:
                data.pop(item)

        project = Dict.from_dict(data, delim='.')
        if project is None:
            project = ws.new_project()
        return project
