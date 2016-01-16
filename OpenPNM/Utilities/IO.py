import scipy as _sp
import numpy as _np
import pandas as _pd
import yaml as _yaml
import itertools as _itertools
from xml.etree import ElementTree as _ET
import OpenPNM
from OpenPNM.Utilities import misc as _misc
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
ctrl = OpenPNM.Base.Controller()


class VTK():
    r"""
    Class for writing a Vtp file to be read by ParaView

    """

    _TEMPLATE = '''
    <?xml version="1.0" ?>
    <VTKFile byte_order="LittleEndian" type="PolyData" version="0.1">
        <PolyData>
            <Piece NumberOfLines="0" NumberOfPoints="0">
                <Points>
                </Points>
                <Lines>
                </Lines>
                <PointData>
                </PointData>
                <CellData>
                </CellData>
            </Piece>
        </PolyData>
    </VTKFile>
    '''.strip()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @staticmethod
    def save(network, filename='', phases=[]):
        r"""
        Save network and phase data to a single vtp file for visualizing in
        Paraview

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network containing the data to be written

        filename : string, optional
            Filename to write data.  If no name is given the file is named
            after ther network

        phases : list, optional
            A list contain OpenPNM Phase object(s) containing data to be
            written

        """

        if filename == '':
            filename = network.name
        filename = filename.split('.')[0] + '.vtp'

        root = _ET.fromstring(VTK._TEMPLATE)
        objs = []
        if type(phases) != list:
            phases = [phases]
        for phase in phases:
            objs.append(phase)
        objs.append(network)
        am = _misc.amalgamate_data(objs=objs)
        key_list = list(sorted(am.keys()))
        points = network['pore.coords']
        pairs = network['throat.conns']

        num_points = len(points)
        num_throats = len(pairs)

        piece_node = root.find('PolyData').find('Piece')
        piece_node.set("NumberOfPoints", str(num_points))
        piece_node.set("NumberOfLines", str(num_throats))

        points_node = piece_node.find('Points')
        coords = VTK._array_to_element("coords", points.T.ravel('F'), n=3)
        points_node.append(coords)

        lines_node = piece_node.find('Lines')
        connectivity = VTK._array_to_element("connectivity", pairs)
        lines_node.append(connectivity)
        offsets = VTK._array_to_element("offsets", 2*_np.arange(len(pairs))+2)
        lines_node.append(offsets)

        point_data_node = piece_node.find('PointData')
        for key in key_list:
            array = am[key]
            if array.dtype == _np.bool:
                array = array.astype(int)
            if array.size != num_points:
                continue
            element = VTK._array_to_element(key, array)
            point_data_node.append(element)

        cell_data_node = piece_node.find('CellData')
        for key in key_list:
            array = am[key]
            if array.dtype == _np.bool:
                array = array.astype(int)
            if array.size != num_throats:
                continue
            element = VTK._array_to_element(key, array)
            cell_data_node.append(element)

        tree = _ET.ElementTree(root)
        tree.write(filename)

        # Make pretty
        with open(filename, 'r+') as f:
            string = f.read()
            string = string.replace('</DataArray>', '</DataArray>\n\t\t\t')
            f.seek(0)
            # consider adding header: '<?xml version="1.0"?>\n'+
            f.write(string)

    @staticmethod
    def load(filename, network=None):
        r"""
        Read in pore and throat data from a saved VTK file.

        Parameters
        ----------
        filename : string
            The name of the 'vtk' file to open.

        Notes
        -----
        This will NOT reproduce original simulation, since all models and
        object relationships are lost.
        """
        return_flag = False
        if network is None:
            network = OpenPNM.Network.Import()
            return_flag = True

        filename = filename.rsplit('.', maxsplit=1)[0]
        tree = _ET.parse(filename+'.vtp')
        piece_node = tree.find('PolyData').find('Piece')

        # Extract connectivity
        conn_element = piece_node.find('Lines').find('DataArray')
        array = VTK._element_to_array(conn_element, 2)
        network['throat.conns'] = array.T

        for element in piece_node.find('PointData').iter('DataArray'):
            key = element.get('Name')
            array = VTK._element_to_array(element)
            netname = key.split('.')[0]
            propname = key.strip(netname+'.')
            network[propname] = array

        if return_flag:
            return network

    @staticmethod
    def _array_to_element(name, array, n=1):
        dtype_map = {
            'int8': 'Int8',
            'int16': 'Int16',
            'int32': 'Int32',
            'int64': 'Int64',
            'uint8': 'UInt8',
            'uint16': 'UInt16',
            'uint32': 'UInt32',
            'uint64': 'UInt64',
            'float32': 'Float32',
            'float64': 'Float64',
            'str': 'String',
        }
        element = _ET.Element('DataArray')
        element.set("Name", name)
        element.set("NumberOfComponents", str(n))
        element.set("type", dtype_map[str(array.dtype)])
        element.text = '\t'.join(map(str, array.ravel()))
        return element

    @staticmethod
    def _element_to_array(element, n=1):
        string = element.text
        dtype = element.get("type")
        array = _np.fromstring(string, sep='\t')
        array = array.astype(dtype)
        if n is not 1:
            array = array.reshape(array.size//n, n)
        return array


class MAT():
    r"""
    Class for reading and writing OpenPNM data to a Matlab 'mat' file
    """

    @staticmethod
    def save(network, filename='', phases=[]):
        r"""
        Write Network to a Mat file for exporting to Matlab. This method will
        be enhanced in a future update, and it's functionality may change!

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
        filename = filename.split('.')[0] + '.mat'
        if phases:  # Ensure it's a list
            phases = list(phases)

        pnMatlab = {}
        new = []
        old = []
        for keys in list(network.keys()):
            old.append(keys)
            new.append(keys.replace('.', '_'))

        for i in range(len(network)):
            pnMatlab[new[i]] = network[old[i]]

        if type(phases) != list:
            phases = [phases]
        if len(phases) != 0:
            for j in range(len(phases)):
                new = []
                old = []

                for keys in list(phases[j].keys()):
                    old.append(keys)
                    new.append(phases[j].name+'_'+keys.replace('.', '_'))

                for i in range(len(phases[j])):
                    pnMatlab[new[i]] = phases[j][old[i]]

        _sp.io.savemat(file_name=filename, mdict=pnMatlab)

    @staticmethod
    def load(filename, network={}, overwrite=True):
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

        overwrite : bool (default is True)
            Indicates whether existing data should be over written if a
            conflicting entry exists in the CSV file.

        Returns
        -------
        If no Network object is supplied then one will be created and returned.

        Notes
        -----
        The 'mat' file must contain data formatted as follows:

        1. The file can contain either or both pore and throat data.

        2. The property names should be in the format of ``pore_volume`` or
        ``throat_surface_area`. In OpenPNM the first \'_\' will be replaced by
        a \'.\' to give \'pore.volume\' or \'throat.surface_area\'.

        3. If pore data is included in the file, then ``pore_coords`` must be
        present, and if throat data is present then ``throat_conns`` must be
        present.

        4. Boolean data represented as 1's and 0's will be converted to the
        Python boolean True and False.  These will become \'labels\' in
        OpenPNM.


        """
        if network == {}:
            network = OpenPNM.Network.Import()
        net = {}

        import scipy.io as _spio
        data = _spio.loadmat(filename)
        # Deal with pore coords and throat conns specially
        if 'throat_conns' in data.keys():
            net.update({'throat.conns': _sp.vstack(data['throat_conns'])})
            Nt = _sp.shape(net['throat.conns'])[0]
            net.update({'throat.all': _sp.ones((Nt,), dtype=bool)})
            del data['throat_conns']
        else:
            logger.warning('\'throat_conns\' not found')
        if 'pore_coords' in data.keys():
            net.update({'pore.coords': _sp.vstack(data['pore_coords'])})
            Np = _sp.shape(net['pore.coords'])[0]
            net.update({'pore.all': _sp.ones((Np,), dtype=bool)})
            del data['pore_coords']
        else:
            logger.warning('\'pore_coords\' not found')

        # Now parse through all the other items
        items = [i for i in data.keys() if '__' not in i]
        for item in items:
            element = item.split('_')[0]
            prop = item.split('_', maxsplit=1)[1]
            vals = _sp.squeeze(data[item].T)
            # If data is not standard array, convert vals appropriately
            if (_sp.sum(vals == 1) + _sp.sum(vals == 0)) \
                    == _sp.shape(vals)[0]:  # If boolean
                vals = vals.astype(bool)
            else:  # If data is an array of lists
                pass
            net[element+'.'+prop] = vals

        network = _update_network(network=network,
                                  net=net,
                                  overwrite=overwrite)
        return network


class Pandas():

    @staticmethod
    def get_data_frames(network, phases=[]):
        r"""
        Convert the Network (and optionally Phase) data to Pandas DataFrames.

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network containing the data to be stored

        phases : list of OpenPNM Phase Objects
            The data on each supplied phase will be added to the CSV file

        Returns
        -------
        A dict containing 2 Pandas DataFrames with 'pore' and 'throat' data in
        each.
        """
        if phases:  # Ensure it's a list
            phases = list(phases)

        # Initialize pore and throat data dictionary with conns and coords
        pdata = {}
        tdata = {}

        # Gather list of prop names from network and geometries
        pprops = set(network.props('pore') + network.labels('pore'))
        for item in network._geometries:
            pprops = pprops.union(set(item.props('pore')))
        tprops = set(network.props('throat') + network.labels('throat'))
        for item in network._geometries:
            tprops = tprops.union(set(item.props('throat')))

        # Select data from network and geometries using keys
        for item in pprops:
            pdata.update({item: network[item]})
        for item in tprops:
            tdata.update({item: network[item]})

        # Gather list of prop names from phases and physics
        for phase in phases:
            # Gather list of prop names
            pprops = set(phase.props('pore') + phase.labels('pore'))
            for item in phase._physics:
                pprops = pprops.union(set(item.props('pore')))
            tprops = set(phase.props('throat') + phase.labels('throat'))
            for item in phase._physics:
                tprops = tprops.union(set(item.props('throat')))
            # Add props to tdata and pdata
            for item in pprops:
                pdata.update({item+'|'+phase.name: phase[item]})
            for item in tprops:
                tdata.update({item+'|'+phase.name: phase[item]})

        # Scan data and convert non-1d arrays to strings
        for item in list(pdata.keys()):
            if _sp.shape(pdata[item]) != (network.Np,):
                array = pdata.pop(item)
                temp = _sp.empty((_sp.shape(array)[0], ), dtype=object)
                for row in range(temp.shape[0]):
                    temp[row] = str(array[row, :]).strip('[]')
                pdata.update({item: temp})

        for item in list(tdata.keys()):
            if _sp.shape(tdata[item]) != (network.Nt,):
                array = tdata.pop(item)
                temp = _sp.empty((_sp.shape(array)[0], ), dtype=object)
                for row in range(temp.shape[0]):
                    temp[row] = str(array[row, :]).strip('[]')
                tdata.update({item: temp})

        data = {'pore.DataFrame': _pd.DataFrame.from_dict(pdata),
                'throat.DataFrame': _pd.DataFrame.from_dict(tdata)}

        return data


class CSV():

    @staticmethod
    def save(network, filename='', phases=[]):
        r"""
        Save all the pore and throat property data on the Network (and
        optionally on any Phases objects) to CSV files.

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network containing the data to be stored

        filename : string
            The name of the file to store the data

        phases : list of OpenPNM Phase Objects
            The data on each supplied phase will be added to the CSV file.

        Notes
        -----
        The data from all Geometry objects is added to the file automatically.
        Furthermore, the Physics data is added for each Phase object that is
        provided.
        """
        if filename == '':
            filename = network.name
        filename = filename.rstrip('.csv')
        if phases:  # Ensure it's a list
            phases = list(phases)
        dataframes = Pandas.get_data_frames(network=network, phases=phases)
        dfp = dataframes['pore.DataFrame']
        dft = dataframes['throat.DataFrame']
        b = dft.join(other=dfp, how='left')
        with open(filename+'.csv', mode='x') as f:
            b.to_csv(f, index=False)

    @staticmethod
    def load(filename, network={}, overwrite=True):
        r"""
        Accepts a file name, reads in the data, and adds it to the Network

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        overwrite : bool (default is True)
            Indicates whether existing data should be over written if a
            conflicting entry exists in the CSV file.

        Notes
        -----
        There are a few rules governing how the data should be stored:

        1. The first row of the file (column headers) must contain the
        property names. The subsequent rows contain the data.

        2. The property names should be in the usual OpenPNM format, such as
        of *pore.volume* or *throat.surface_area*.

        3. Each column represents a specific property.  For Np x 1 or Nt x 1
        data such as *pore.volume* this is straightforward.  For Np x m or
        Nt x m data, it must be entered in as a set of values NOT separated by
        commas.  For instance, the *pore.coords* values should be X Y Z with
        spaces, not commas between them.

        4. The file can contain both or either pore and throat data.

        5. Labels can be imported by placing the characters TRUE and FALSE
        in a column corresponding to the label name (i.e. *pore.front*).  TRUE
        indicates where the label applies and FALSE otherwise.
        """
        if network == {}:
            network = OpenPNM.Network.Import()
        # Instantiate new empty dict
        net = {}

        filename = filename.rstrip('.csv') + '.csv'
        with open(filename) as f:
            a = _pd.read_table(filepath_or_buffer=f,
                               sep=',',
                               skipinitialspace=True,
                               index_col=False,
                               true_values=['T', 't', 'True', 'true',
                                            'TRUE'],
                               false_values=['F', 'f', 'False', 'false',
                                             'FALSE'])

        # Now parse through all the other items
        for item in a.keys():
            element = item.split('.')[0]
            prop = item.split('.', maxsplit=1)[1]
            data = _sp.array(a[item].dropna())
            if type(data[0]) is str:
                N = _sp.shape(data)[0]
                if '.' in data[0].split(' ')[0]:  # Decimal means float
                    dtype = float
                else:
                    dtype = int
                temp = _sp.empty(_sp.shape(data), dtype=object)
                for row in range(N):
                    temp[row] = _sp.fromstring(data[row], sep=' ', dtype=dtype)
                data = _sp.vstack(temp)
            else:
                dtype = type(data[0])
            net[element+'.'+prop] = data.astype(dtype)

        network = _update_network(network=network,
                                  net=net,
                                  overwrite=overwrite)
        return network


class YAML():
    r"""
    This format is meant specifcally for exchanging data with NetworkX, which
    is a common tool for dealing with network structures.
    """

    @staticmethod
    def save():
        # TODO: This would be a great place for a new developer to contribute
        raise NotImplemented

    @staticmethod
    def load(filename, network={}, overwrite=True):
        r"""
        Add data to an OpenPNM Network from a NetworkX generated YAML file.

        Parameters
        ----------
        filename : string
            The yaml file containing the NetworkX data

        network : OpenPNM Network Object
            The OpenPNM Network onto which the data should be loaded.  If no
            Network is supplied then an empty Import Network is created and
            returned.

        overwite : boolean
            A flag to indicate whether data existing on the Network should be
            overwritten by data in the file or not.  This is useful when adding
            data in a sequential manner from several different sources.

        """
        if network == {}:
            network = OpenPNM.Network.Import()
        # Instantiate new empty dict
        net = {}

        # Open file and read first line, to prevent networkx instantiation
        with open(filename) as f:
            line = f.readline()
            if line.starts('!!python/object:networkx.classes.graph.Graph'):
                a = _yaml.safe_load(f)
            else:
                raise ('Provided file does not appear to a NetworkX file')

        # Parsing node data
        Np = len(a['node'])
        net.update({'pore.all': _sp.ones((Np,), dtype=bool)})
        for n in a['node'].keys():
            props = a['node'][n]
            for item in props.keys():
                val = a['node'][n][item]
                if 'pore.'+item not in net.keys():
                    dtype = type(val)
                    if dtype is list:
                        dtype = type(val[0])
                        cols = len(val)
                        net['pore.'+item] = _sp.ndarray((Np, cols),
                                                        dtype=dtype)
                    else:
                        net['pore.'+item] = _sp.ndarray((Np,), dtype=dtype)
                net['pore.'+item][n] = val

        # Parsing edge data
        # Deal with conns explicitly
        conns = []
        for n in a['edge'].keys():
            neighbors = a['edge'][n].keys()
            conns.extend([sorted([i, n]) for i in neighbors])
        # Remove duplicate pairs from conns and sort
        conns.sort()
        conns = list(conns for conns, _ in _itertools.groupby(conns))
        # Add conns to Network
        Nt = len(conns)
        net.update({'throat.all': _sp.ones(Nt, dtype=bool)})
        net.update({'throat.conns': _sp.array(conns)})

        # Scan through each edge and extract all its properties
        i = 0
        for t in conns:
            props = a['edge'][t[0]][t[1]]
            for item in props:
                val = props[item]
                if 'throat.'+item not in net.keys():
                    dtype = type(val)
                    if dtype is list:
                        dtype = type(val[0])
                        cols = len(val)
                        net['throat.'+item] = _sp.ndarray((Nt, cols),
                                                          dtype=dtype)
                    else:
                        net['throat.'+item] = _sp.ndarray((Nt,),
                                                          dtype=dtype)
                net['throat.'+item][i] = val
            i += 1

        network = _update_network(network=network,
                                  net=net,
                                  overwrite=overwrite)
        return network


def _update_network(network, net, overwrite):
    # Add newly read props to the network
    Np = [_sp.shape(net[i])[0] for i in net.keys() if i.startswith('pore')]
    if Np and _sp.all(Np == Np[0]):
        net.update({'pore.all': _sp.ones((Np[0],), dtype=bool)})
    Nt = [_sp.shape(net[i])[0] for i in net.keys() if i.startswith('throat')]
    if Nt and _sp.all(Nt == Nt[0]):
        net.update({'throat.all': _sp.ones((Nt[0],), dtype=bool)})
    for item in net.keys():
        if overwrite:
            network.update({item: net[item]})
        elif item not in network:
            network.update({item: net[item]})
        else:
            logger.warning('\''+item+'\' already present')
