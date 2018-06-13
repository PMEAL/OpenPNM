import os as os
import scipy as sp
import pandas as pd
from openpnm.topotools import trim
from openpnm.core import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
from pathlib import Path
logger = logging.getLogger(__name__)


class Statoil(GenericIO):
    r"""
    This class is for loading data stored in the 'Statoil' file format.  More
    specifically, this file format is used by the network extraction code of
    Blunt's group at Imperial College London, so this class can be used to load
    and work with those network.  Numerous datasets are available for download
    from the group's `website <http://tinyurl.com/zurko4q>`_.

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
        net = {}

        # ---------------------------------------------------------------------
        # Parse the link1 file
        path = Path(path)
        filename = Path(path.resolve(), prefix+'_link1.dat')
        with open(filename, mode='r') as f:
            link1 = pd.read_table(filepath_or_buffer=f,
                                  header=None,
                                  skiprows=1,
                                  sep=' ',
                                  skipinitialspace=True,
                                  index_col=0)
        link1.columns = ['throat.pore1', 'throat.pore2', 'throat.radius',
                         'throat.shape_factor', 'throat.total_length']
        # Add link1 props to net
        net['throat.conns'] = sp.vstack((link1['throat.pore1']-1,
                                         link1['throat.pore2']-1)).T
        net['throat.conns'] = sp.sort(net['throat.conns'], axis=1)
        net['throat.radius'] = sp.array(link1['throat.radius'])
        net['throat.shape_factor'] = sp.array(link1['throat.shape_factor'])
        net['throat.total_length'] = sp.array(link1['throat.total_length'])
        # ---------------------------------------------------------------------
        filename = Path(path.resolve(), prefix+'_link2.dat')
        with open(filename, mode='r') as f:
            link2 = pd.read_table(filepath_or_buffer=f,
                                  header=None,
                                  sep=' ',
                                  skipinitialspace=True,
                                  index_col=0)
        link2.columns = ['throat.pore1', 'throat.pore2',
                         'throat.pore1_length', 'throat.pore2_length',
                         'throat.length', 'throat.volume',
                         'throat.clay_volume']
        # Add link2 props to net
        net['throat.length'] = sp.array(link2['throat.length'])
        net['throat.volume'] = sp.array(link2['throat.volume'])
        net['throat.clay_volume'] = sp.array(link2['throat.clay_volume'])
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
        node1 = pd.DataFrame(array[:, [1, 2, 3, 4]])
        node1.columns = ['pore.x_coord', 'pore.y_coord', 'pore.z_coord',
                         'pore.coordination_number']
        # Add node1 props to net
        net['pore.coords'] = sp.vstack((node1['pore.x_coord'],
                                        node1['pore.y_coord'],
                                        node1['pore.z_coord'])).T
        # ---------------------------------------------------------------------
        # Parse the node1 file
        filename = Path(path.resolve(), prefix+'_node2.dat')
        with open(filename, mode='r') as f:
            node2 = pd.read_table(filepath_or_buffer=f,
                                  header=None,
                                  sep=' ',
                                  skipinitialspace=True,
                                  index_col=0)
        node2.columns = ['pore.volume', 'pore.radius', 'pore.shape_factor',
                         'pore.clay_volume']
        # Add node2 props to net
        net['pore.volume'] = sp.array(node2['pore.volume'])
        net['pore.radius'] = sp.array(node2['pore.radius'])
        net['pore.shape_factor'] = sp.array(node2['pore.shape_factor'])
        net['pore.clay_volume'] = sp.array(node2['pore.clay_volume'])

        if network is None:
            network = GenericNetwork()
        network = cls._update_network(network=network, net=net)

        # Use OpenPNM Tools to clean up network
        # Trim throats connected to 'inlet' or 'outlet' reservoirs
        trim1 = sp.where(sp.any(net['throat.conns'] == -1, axis=1))[0]
        # Apply 'outlet' label to these pores
        outlets = network['throat.conns'][trim1, 1]
        network['pore.outlets'] = False
        network['pore.outlets'][outlets] = True
        trim2 = sp.where(sp.any(net['throat.conns'] == -2, axis=1))[0]
        # Apply 'inlet' label to these pores
        inlets = network['throat.conns'][trim2, 1]
        network['pore.inlets'] = False
        network['pore.inlets'][inlets] = True
        # Now trim the throats
        to_trim = sp.hstack([trim1, trim2])
        trim(network=network, throats=to_trim)

        return network.project
