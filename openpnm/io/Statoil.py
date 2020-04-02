import numpy as np
import scipy as sp
from openpnm.topotools import trim
from openpnm.utils import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
from pathlib import Path
logger = logging.getLogger(__name__)


class Statoil(GenericIO):
    r"""
    The StatOil format is used by the Maximal Ball network extraction code of
    the Imperial College London group

    This class can be used to load and work with those networks.  Numerous
    datasets are available for download from the group's
    `website <http://tinyurl.com/zurko4q>`_.

    The so-called 'Statoil' format consists of 4 different files in a single
    folder.  The data is stored in columns with each corresponding to a
    specific property.  Headers are not provided in the files, so one must
    refer to various theses and documents to interpret their meaning.
    """
    @classmethod
    def load(cls, path, prefix, network=None):
        r"""
        Load data from the \'dat\' files located in specified folder.

        Parameters
        ----------
        path : string
            The full path to the folder containing the set of \'dat\' files.

        prefix : string
            The file name prefix on each file. The data files are stored
            as \<prefix\>_node1.dat.

        network : OpenPNM Network Object
            If given then the data will be loaded on it and returned.  If not
            given, a Network will be created and returned.

        Returns
        -------
        An OpenPNM Project containing a GenericNetwork holding all the data

        """
        from pandas import read_table, DataFrame

        net = {}

        # Parse the link1 file
        path = Path(path)
        filename = Path(path.resolve(), prefix+'_link1.dat')
        with open(filename, mode='r') as f:
            link1 = read_table(filepath_or_buffer=f,
                                  header=None,
                                  skiprows=1,
                                  sep=' ',
                                  skipinitialspace=True,
                                  index_col=0)
        link1.columns = ['throat.pore1', 'throat.pore2', 'throat.radius',
                         'throat.shape_factor', 'throat.total_length']
        # Add link1 props to net
        net['throat.conns'] = np.vstack((link1['throat.pore1']-1,
                                         link1['throat.pore2']-1)).T
        net['throat.conns'] = np.sort(net['throat.conns'], axis=1)
        net['throat.radius'] = np.array(link1['throat.radius'])
        net['throat.shape_factor'] = np.array(link1['throat.shape_factor'])
        net['throat.total_length'] = np.array(link1['throat.total_length'])

        filename = Path(path.resolve(), prefix+'_link2.dat')
        with open(filename, mode='r') as f:
            link2 = read_table(filepath_or_buffer=f,
                                  header=None,
                                  sep=' ',
                                  skipinitialspace=True,
                                  index_col=0)
        link2.columns = ['throat.pore1', 'throat.pore2',
                         'throat.pore1_length', 'throat.pore2_length',
                         'throat.length', 'throat.volume',
                         'throat.clay_volume']
        # Add link2 props to net
        cl_t = np.array(link2['throat.length'])
        net['throat.length'] = cl_t
        net['throat.conduit_lengths.throat'] = cl_t
        net['throat.volume'] = np.array(link2['throat.volume'])
        cl_p1 = np.array(link2['throat.pore1_length'])
        net['throat.conduit_lengths.pore1'] = cl_p1
        cl_p2 = np.array(link2['throat.pore2_length'])
        net['throat.conduit_lengths.pore2'] = cl_p2
        net['throat.clay_volume'] = np.array(link2['throat.clay_volume'])
        # ---------------------------------------------------------------------
        # Parse the node1 file
        filename = Path(path.resolve(), prefix+'_node1.dat')
        with open(filename, mode='r') as f:
            row_0 = f.readline().split()
            num_lines = int(row_0[0])
            array = sp.ndarray([num_lines, 6])
            for i in range(num_lines):
                row = f.readline()\
                       .replace('\t', ' ').replace('\n', ' ').split()
                array[i, :] = row[0:6]
        node1 = DataFrame(array[:, [1, 2, 3, 4]])
        node1.columns = ['pore.x_coord', 'pore.y_coord', 'pore.z_coord',
                         'pore.coordination_number']
        # Add node1 props to net
        net['pore.coords'] = np.vstack((node1['pore.x_coord'],
                                        node1['pore.y_coord'],
                                        node1['pore.z_coord'])).T
        # ---------------------------------------------------------------------
        # Parse the node1 file
        filename = Path(path.resolve(), prefix+'_node2.dat')
        with open(filename, mode='r') as f:
            node2 = read_table(filepath_or_buffer=f,
                                  header=None,
                                  sep=' ',
                                  skipinitialspace=True,
                                  index_col=0)
        node2.columns = ['pore.volume', 'pore.radius', 'pore.shape_factor',
                         'pore.clay_volume']
        # Add node2 props to net
        net['pore.volume'] = np.array(node2['pore.volume'])
        net['pore.radius'] = np.array(node2['pore.radius'])
        net['pore.shape_factor'] = np.array(node2['pore.shape_factor'])
        net['pore.clay_volume'] = np.array(node2['pore.clay_volume'])
        net['throat.area'] = ((net['throat.radius']**2)
                              / (4.0*net['throat.shape_factor']))
        net['pore.area'] = ((net['pore.radius']**2)
                            / (4.0*net['pore.shape_factor']))

        if network is None:
            network = GenericNetwork()
        network = cls._update_network(network=network, net=net)

        # Use OpenPNM Tools to clean up network
        # Trim throats connected to 'inlet' or 'outlet' reservoirs
        trim1 = np.where(np.any(net['throat.conns'] == -1, axis=1))[0]
        # Apply 'outlet' label to these pores
        outlets = network['throat.conns'][trim1, 1]
        network['pore.outlets'] = False
        network['pore.outlets'][outlets] = True
        trim2 = np.where(np.any(net['throat.conns'] == -2, axis=1))[0]
        # Apply 'inlet' label to these pores
        inlets = network['throat.conns'][trim2, 1]
        network['pore.inlets'] = False
        network['pore.inlets'][inlets] = True
        # Now trim the throats
        to_trim = np.hstack([trim1, trim2])
        trim(network=network, throats=to_trim)

        return network.project
