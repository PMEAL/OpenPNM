import os
import numpy as np
from openpnm.topotools import trim, extend
from openpnm.utils import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
from openpnm.geometry import GenericGeometry
import openpnm.models as mods
from pathlib import Path
from pandas import read_table, DataFrame
from tqdm import tqdm
logger = logging.getLogger(__name__)


class Statoil(GenericIO):
    r"""
    The StatOil format is used by the Maximal Ball network extraction code of
    the Imperial College London group

    This class can be used to load and work with those networks. Numerous
    datasets are available for download from the group's
    `website <http://tinyurl.com/zurko4q>`_.

    The 'Statoil' format consists of 4 different files in a single
    folder. The data is stored in columns with each corresponding to a
    specific property. Headers are not provided in the files, so one must
    refer to various theses and documents to interpret their meaning.

    """

    @classmethod
    def export_data(cls, network, shape, prefix=None, path=None, Pin=None, Pout=None):
        r"""

        Parameters
        ----------
        network : GenericNetwork
            The network
        shape : array_like
            An ndim-by-1 array or list containing the network dimensions
            in physical units (i.e. um)
        prefix : str
            The prefix to append to each file name, such as
            ``<prefix>_node1.dat``. If not provided ``project.name`` is used.
        path : str or path object
            The location where the exported files should be stored. If not
            provided the current working directory is used
        Pinlet and Poutlet : scalar, int (optional)
            The pore index of the inlet and outlet reservoir pores. If not
            provided then it is assumed they are the second last and last
            pores in the network, respectively.  This would be the case if
            the ``add_reservoir_pore`` function had been called prior to
            exporting.

        """
        if path is None:
            path = os.getcwd()
        p = Path(path)
        if prefix is None:
            prefix = network.project.name
        # Deal with reservoir pores
        if Pin is None:
            Pin = network.Np - 2
        if Pout is None:
            Pout = network.Np - 1
        Ptemp = network.find_neighbor_pores(pores=Pin)
        inlets = np.zeros_like(network.Ps, dtype=bool)
        inlets[Ptemp] = True
        Ptemp = network.find_neighbor_pores(pores=Pout)
        outlets = np.zeros_like(network.Ps, dtype=bool)
        outlets[Ptemp] = True

        # Write link 1 file
        props = ['throat.conns',
                 'throat.diameter',
                 'throat.shape_factor',
                 'throat.total_length']
        with open(p.joinpath(prefix + '_link1.dat'), 'wt') as f:
            f.write(str(network.Nt) + '\n')
            for row in tqdm(network.throats(), desc='Writing Link1 file'):
                s = ''
                s = s + '{:>9}'.format(str(row+1))
                for col in props:
                    try:
                        val = network[col][row]
                    except KeyError:
                        val = 0.0
                    if col == 'throat.conns':
                        val = np.copy(val)
                        val[val == network.Np - 1] = -1
                        val[val == (network.Np - 2)] = -2
                        s = s + '{:>9}'.format(str(val[0] + 1))
                        s = s + '{:>9}'.format(str(val[1] + 1))
                        continue
                    if isinstance(val, float):
                        if 'diameter' in col:  # Convert to radius
                            val = val/2
                        if np.isnan(val):
                            val = 0.0
                        val = np.format_float_scientific(val, precision=6,
                                                         exp_digits=3,
                                                         trim='k',
                                                         unique=False)
                    s = s + '{:>15}'.format(str(val))
                s = s + '\n'  # Remove trailing tab and add a new line
                f.write(s)

        # Write Link 2 file
        props = ['throat.conns',
                 'throat.conduit_lengths.pore1',
                 'throat.conduit_lengths.pore2',
                 'throat.conduit_lengths.throat',
                 'throat.volume',
                 'throat.clay_volume']
        with open(p.joinpath(prefix + '_link2.dat'), 'wt') as f:
            for row in tqdm(network.throats(), desc='Writing Link2 file'):
                s = ''
                s = s + '{:>9}'.format(str(row+1))
                for col in props:
                    try:
                        val = network[col][row]
                    except KeyError:
                        val = 0.0
                    if col == 'throat.conns':
                        val = np.copy(val)
                        val[val == network.Np - 1] = -1
                        val[val == (network.Np - 2)] = -2
                        # Original file has 7 spaces for pore indices, but
                        # this is not enough for networks with > 10 million
                        # pores so I have bumped it to 9. I'm not sure if
                        # this will still work with the ICL binaries.
                        s = s + '{:>9}'.format(str(val[0] + 1))
                        s = s + '{:>9}'.format(str(val[1] + 1))
                        continue
                    if isinstance(val, float):
                        if np.isnan(val):
                            val = 0.0
                        val = np.format_float_scientific(val, precision=6,
                                                         exp_digits=3,
                                                         trim='k',
                                                         unique=False)
                    s = s + '{:>15}'.format(str(val))
                s = s + '\n'  # Remove trailing tab and a new line
                f.write(s)

        # Write Node 1 file
        with open(p.joinpath(prefix + '_node1.dat'), 'wt') as f:
            s = ''
            s = s + str(network.num_pores('reservoir', 'not'))
            for d in shape:
                val = np.format_float_scientific(d, precision=6, exp_digits=3,
                                                 trim='k', unique=False)
                s = s + '{:>17}'.format(str(val))
            s = s + '\n'
            f.write(s)
            for row in tqdm(network.pores('reservoir', mode='not'),
                            desc='Writing Node1 file'):
                if row in [Pin, Pout]:
                    continue
                s = ''
                s = s + '{:>9}'.format(str(row+1))
                for c in network['pore.coords'][row]:
                    if isinstance(c, float):
                        c = np.format_float_scientific(c, precision=6,
                                                       exp_digits=3,
                                                       trim='k',
                                                       unique=False)
                    s = s + '{:>15}'.format(str(c))
                s = s + '{:>9}'.format(str(network.num_neighbors(row)[0]))
                for n in network.find_neighbor_pores(row):
                    if n == Pin:
                        n = 0
                    elif n == Pout:
                        n = -1
                    else:
                        n = n + 1
                    s = s + '{:>9}'.format(str(n))
                s = s + '{:>9}'.format(str(int(inlets[row])))
                s = s + '{:>9}'.format(str(int(outlets[row])))
                for n in network.find_neighbor_throats(row):
                    s = s + '{:>9}'.format(str(n + 1))
                s = s + '\n'  # Remove trailing tab and a new line
                f.write(s)

        # Write Node 2 file
        props = ['pore.volume',
                 'pore.diameter',
                 'pore.shape_factor',
                 'pore.clay_volume']
        with open(p.joinpath(prefix + '_node2.dat'), 'wt') as f:
            for row in tqdm(network.pores('reservoir', mode='not'),
                            desc='Writing Node2 file'):
                s = ''
                s = s + '{:>9}'.format(str(row+1))
                for col in props:
                    try:
                        val = network[col][row]
                    except KeyError:
                        val = 0.0
                    if isinstance(val, float):
                        if 'diameter' in col:
                            val = val/2
                        val = np.format_float_scientific(val, precision=6,
                                                         exp_digits=3,
                                                         trim='k',
                                                         unique=False)
                        s = s + '{:>15}'.format(str(val))
                s = s + '\n'  # Remove trailing tab and a new line
                f.write(s)


    @classmethod
    def import_data(cls, path, prefix, network=None):
        r"""
        Load data from the \'dat\' files located in specified folder.

        Parameters
        ----------
        path : str
            The full path to the folder containing the set of \'dat\' files.
        prefix : str
            The file name prefix on each file. The data files are stored
            as \<prefix\>_node1.dat.
        network : GenericNetwork
            If given then the data will be loaded on it and returned.  If not
            given, a Network will be created and returned.

        Returns
        -------
        Project
            An OpenPNM Project containing a GenericNetwork holding all the data

        """
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
            array = np.ndarray([num_lines, 6])
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
        # Parse the node2 file
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


def from_statoil(path, prefix, network=None):
    Statoil.import_data(path=path, prefix=prefix, network=network)


from_statoil.__doc__ = Statoil.import_data.__doc__


def to_statoil(network, shape, prefix=None, path=None, Pin=None, Pout=None):
    project = Statoil.export_data(network=network, shape=shape,
                                  prefix=prefix, path=path,
                                  Pin=Pin, Pout=Pout)
    return project


to_statoil.__doc__ = Statoil.import_data.__doc__
