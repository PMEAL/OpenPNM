import scipy as sp
from openpnm.core import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
logger = logging.getLogger(__name__)


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
    def save(cls, network, phases=[], filename=''):
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
        if filename == '':
            filename = network.name
        filename = filename.replace('.mat', '') + '.mat'
        if type(phases) is not list:  # Ensure it's a list
            phases = [phases]

        keys = network.props(deep=True) + network.labels()
        pnMatlab = {i.replace('.', '_'): network[i] for i in keys}

        for phase in phases:
            keys = phase.props(mode=['all', 'deep']) + phase.labels()
            temp = {i.replace('.', '_')+'|'+phase.name: phase[i]
                    for i in keys}
            pnMatlab.update(temp)

        sp.io.savemat(file_name=filename, mdict=pnMatlab)

    @classmethod
    def load(cls, filename, network=None, return_geometry=False):
        r"""
        Loads data onto the given network from an appropriately formatted
        'mat' file (i.e. MatLAB output).

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        network : OpenPNM Network Object
            The Network object onto which the data should be loaded.  If no
            Network is supplied than one will be created and returned.

        return_geometry : Boolean
            If True, then all geometrical related properties are removed from
            the Network object and added to a GenericGeometry object.  In this
            case the method returns a tuple containing (network, geometry). If
            False (default) then the returned Network will contain all
            properties that were in the original file.  In this case, the user
            can call the ```split_geometry``` method explicitly to perform the
            separation.

        Returns
        -------
        If no Network object is supplied then one will be created and returned.

        If return_geometry is True, then a tuple is returned containing both
        the network and a geometry object.

        """
        net = {}

        import scipy.io as spio
        data = spio.loadmat(filename)
        # Deal with pore coords and throat conns specially
        if 'throat_conns' in data.keys():
            net.update({'throat.conns': sp.vstack(data['throat_conns'])})
            Nt = sp.shape(net['throat.conns'])[0]
            net.update({'throat.all': sp.ones((Nt,), dtype=bool)})
            del data['throat_conns']
        else:
            logger.warning('\'throat_conns\' not found')
        if 'pore_coords' in data.keys():
            net.update({'pore.coords': sp.vstack(data['pore_coords'])})
            Np = sp.shape(net['pore.coords'])[0]
            net.update({'pore.all': sp.ones((Np,), dtype=bool)})
            del data['pore_coords']
        else:
            logger.warning('\'pore_coords\' not found')

        # Now parse through all the other items
        items = [i for i in data.keys() if '__' not in i]
        for item in items:
            element = item.split('_')[0]
            prop = item.split('_', maxsplit=1)[1]
            net[element+'.'+prop] = sp.squeeze(data[item].T)

        if network is None:
            network = GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)
        return network
