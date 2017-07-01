import os as _os
import itertools as _itertools
from xml.etree import ElementTree as _ET
import scipy as _sp
import numpy as _np
import pandas as _pd
import yaml as _yaml
import OpenPNM
from OpenPNM.Utilities import misc as _misc
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
mgr = OpenPNM.Base.Workspace()


class GenericIO():

    @classmethod
    def save(cls):
        raise NotImplementedError("The \'save\' method for this class " +
                                  "does not exist yet")

    @classmethod
    def load(cls):
        raise NotImplementedError("The \'load\' method for this class " +
                                  "does not exist yet")

    @staticmethod
    def split_geometry(network):
        r"""
        This method accepts an OpenPNM Network object and removes all geometry
        related pore and throat properties, (basically all values other than
        ```'pore.coords'``` and ```throat.conns```), and places them on a
        GenericGeometry object.  Any labels on the Network are left intact.

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network that possesses the geometrical values

        Returns
        -------
        geometry : OpenPNM Geometry Object
            The new GenericGeometry object that was created to contain the
            geometrical pore and throat properties.

        """
        geom = OpenPNM.Geometry.GenericGeometry(network=network,
                                                pores=network.Ps,
                                                throats=network.Ts)
        for item in network.props():
            if item not in ['pore.coords', 'throat.conns']:
                geom.update({item: network.pop(item)})
        return geom

    @classmethod
    def _update_network(cls, network, net, return_geometry=False):
        # Infer Np and Nt from length of given prop arrays in file
        for element in ['pore', 'throat']:
            N = [_sp.shape(net[i])[0] for i in net.keys() if i.startswith(element)]
            if N:
                N = _sp.array(N)
                if _sp.all(N == N[0]):
                    if (network._count(element) == N[0]) \
                            or (network._count(element) == 0):
                        network.update({element+'.all': _sp.ones((N[0],),
                                                                 dtype=bool)})
                        net.pop(element+'.all', None)
                    else:
                        raise Exception('Length of '+element+' data in file' +
                                        ' does not match network')
                else:
                    raise Exception(element+' data in file have inconsistent' +
                                    ' lengths')
        # Add data on dummy net to actual network
        for item in net.keys():
            # Try to infer array types and change if necessary
            # Chcek for booleans disguised and 1's and 0's
            num0s = _sp.sum(net[item] == 0)
            num1s = _sp.sum(net[item] == 1)
            if (num1s + num0s) == _sp.shape(net[item])[0]:
                net[item] = net[item].astype(bool)
            # Write data to network object
            if item not in network:
                network.update({item: net[item]})
            else:
                logger.warning('\''+item+'\' already present')
        if return_geometry:
            geometry = cls.split_geometry(network)
            network = (network, geometry)
        return network

    @classmethod
    def _write_file(cls, filename, ext):
        ext = ext.replace('.', '').lower()
        if ext not in ['csv', 'yaml', 'mat', 'vtp', 'dat']:
            raise Exception(ext+' is not a supported file extension')
        filename = filename.rstrip('.'+ext)
        filename = filename+'.'+ext
        try:
            logger.warning(filename+' already exists, contents will be ' +
                           'overwritten')
            f = open(filename, mode='w')
        except:
            f = open(filename, mode='x')
        return f

    @classmethod
    def _read_file(cls, filename, ext):
        ext = ext.replace('.', '').lower()
        if ext not in ['csv', 'yaml', 'mat', 'vtp', 'dat']:
            raise Exception(ext+' is not a supported file extension')
        if not filename.endswith('.'+ext):
            filename = filename+'.'+ext
        f = open(filename, mode='r')
        return f


class VTK(GenericIO):
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

    @classmethod
    def save(cls, network, filename='', phases=[], legacy=True):
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

        legacy : boolean
            If True (default) the property names will be of the format
            \'pore.Cubic_asd43_diameter'\, while if False they will be
            \'pore.diameter|Cubic_asd43\'.  The latter style is consistent
            with all of the other IO methods, while the former is compatible
            with existing code, such as Paraview State files.   Eventually,
            this option will be derprecated and removed.

        """

        if filename == '':
            filename = network.name
        if ~filename.endswith('.vtp'):
            filename = filename+'.vtp'

        root = _ET.fromstring(VTK._TEMPLATE)
        objs = []
        if type(phases) != list:
            phases = [phases]
        for phase in phases:
            objs.append(phase)
        objs.append(network)
        if legacy:
            am = _misc.amalgamate_data(objs=objs)
        else:
            am = {i: network[i] for i in
                  network.props(mode=['all', 'deep']) + network.labels()}
            for phase in phases:
                dict_ = {i+'|'+phase.name: phase[i] for i in
                         phase.props(mode=['all', 'deep']) + phase.labels()}
                am.update(dict_)
        key_list = list(sorted(am.keys()))
        points = network['pore.coords']
        pairs = network['throat.conns']

        num_points = _sp.shape(points)[0]
        num_throats = _sp.shape(pairs)[0]

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

    @classmethod
    def load(cls, filename, network=None, return_geometry=False):
        r"""
        Read in pore and throat data from a saved VTK file.

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

        filename = filename.rsplit('.', maxsplit=1)[0]
        tree = _ET.parse(filename+'.vtp')
        piece_node = tree.find('PolyData').find('Piece')

        # Extract connectivity
        conn_element = piece_node.find('Lines').find('DataArray')
        array = VTK._element_to_array(conn_element, 2)
        net.update({'throat.conns': array})
        # Extract coordinates
        coord_element = piece_node.find('Points').find('DataArray')
        array = VTK._element_to_array(coord_element, 3)
        net.update({'pore.coords': array})

        # Extract pore data
        for item in piece_node.find('PointData').iter('DataArray'):
            key = item.get('Name')
            element = key.split('.')[0]
            array = VTK._element_to_array(item)
            propname = key.split('.')[1]
            net.update({element+'.'+propname: array})
        # Extract throat data
        for item in piece_node.find('CellData').iter('DataArray'):
            key = item.get('Name')
            element = key.split('.')[0]
            array = VTK._element_to_array(item)
            propname = key.split('.')[1]
            net.update({element+'.'+propname: array})

        if network is None:
            network = OpenPNM.Network.GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)
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
    def load(cls, path, prefix, network=None, return_geometry=False):
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
            given, a Network will be created and return.

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

        # ---------------------------------------------------------------------
        # Parse the link1 file
        for item in ['link1']:
            filename = _os.path.join(path, prefix+'_'+item+'.dat')
            with cls._read_file(filename=filename, ext='dat') as f:
                link1 = _pd.read_table(filepath_or_buffer=f,
                                       header=None,
                                       skiprows=1,
                                       sep=' ',
                                       skipinitialspace=True,
                                       index_col=0)
        link1.columns = ['throat.pore1', 'throat.pore2', 'throat.radius',
                         'throat.shape_factor', 'throat.total_length']
        # Add link1 props to net
        net['throat.conns'] = _sp.vstack((link1['throat.pore1']-1,
                                          link1['throat.pore2']-1)).T
        net['throat.conns'] = _sp.sort(net['throat.conns'], axis=1)
        net['throat.radius'] = _sp.array(link1['throat.radius'])
        net['throat.shape_factor'] = _sp.array(link1['throat.shape_factor'])
        net['throat.total_length'] = _sp.array(link1['throat.total_length'])
        # ---------------------------------------------------------------------
        # Parse the link2 file
        for item in ['link2']:
            filename = _os.path.join(path, prefix+'_'+item+'.dat')
            with cls._read_file(filename=filename, ext='dat') as f:
                link2 = _pd.read_table(filepath_or_buffer=f,
                                       header=None,
                                       sep=' ',
                                       skipinitialspace=True,
                                       index_col=0)
        link2.columns = ['throat.pore1', 'throat.pore2',
                         'throat.pore1_length', 'throat.pore2_length',
                         'throat.length', 'throat.volume',
                         'throat.clay_volume']
        # Add link2 props to net
        net['throat.length'] = _sp.array(link2['throat.length'])
        net['throat.volume'] = _sp.array(link2['throat.volume'])
        net['throat.clay_volume'] = _sp.array(link2['throat.clay_volume'])
        # ---------------------------------------------------------------------
        # Parse the node1 file
        for item in ['node1']:
            filename = _os.path.join(path, prefix+'_'+item+'.dat')
            with cls._read_file(filename=filename, ext='dat') as f:
                row_0 = f.readline().split()
                num_lines = int(row_0[0])
                array = _sp.ndarray([num_lines, 6])
                for i in range(num_lines):
                    row = f.readline()\
                           .replace('\t', ' ').replace('\n', ' ').split()
                    array[i, :] = row[0:6]
        node1 = _pd.DataFrame(array[:, [1, 2, 3, 4]])
        node1.columns = ['pore.x_coord', 'pore.y_coord', 'pore.z_coord',
                         'pore.coordination_number']
        # Add node1 props to net
        net['pore.coords'] = _sp.vstack((node1['pore.x_coord'],
                                         node1['pore.y_coord'],
                                         node1['pore.z_coord'])).T
        # ---------------------------------------------------------------------
        # Parse the node1 file
        for item in ['node2']:
            filename = _os.path.join(path, prefix+'_'+item+'.dat')
            with cls._read_file(filename=filename, ext='dat') as f:
                node2 = _pd.read_table(filepath_or_buffer=f,
                                       header=None,
                                       sep=' ',
                                       skipinitialspace=True,
                                       index_col=0)
        node2.columns = ['pore.volume', 'pore.radius', 'pore.shape_factor',
                         'pore.clay_volume']
        # Add node2 props to net
        net['pore.volume'] = _sp.array(node2['pore.volume'])
        net['pore.radius'] = _sp.array(node2['pore.radius'])
        net['pore.shape_factor'] = _sp.array(node2['pore.shape_factor'])
        net['pore.clay_volume'] = _sp.array(node2['pore.clay_volume'])

        if network is None:
            network = OpenPNM.Network.GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)

        # Use OpenPNM Tools to clean up network
        # Trim throats connected to 'inlet' or 'outlet' reservoirs
        trim1 = _sp.where(_sp.any(net['throat.conns'] == -1, axis=1))[0]
        # Apply 'outlet' label to these pores
        outlets = network['throat.conns'][trim1, 1]
        network['pore.outlets'] = False
        network['pore.outlets'][outlets] = True
        trim2 = _sp.where(_sp.any(net['throat.conns'] == -2, axis=1))[0]
        # Apply 'inlet' label to these pores
        inlets = network['throat.conns'][trim2, 1]
        network['pore.inlets'] = False
        network['pore.inlets'][inlets] = True
        # Now trim the throats
        trim = _sp.hstack([trim1, trim2])
        network.trim(throats=trim)

        return network


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
    def save(cls, network, filename='', phases=[]):
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

        keys = network.props(mode=['all', 'deep']) + network.labels()
        pnMatlab = {i.replace('.', '_'): network[i] for i in keys}

        for phase in phases:
            keys = phase.props(mode=['all', 'deep']) + phase.labels()
            temp = {i.replace('.', '_')+'|'+phase.name: phase[i]
                    for i in keys}
            pnMatlab.update(temp)

        _sp.io.savemat(file_name=filename, mdict=pnMatlab)

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
            net[element+'.'+prop] = _sp.squeeze(data[item].T)

        if network is None:
            network = OpenPNM.Network.GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)
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
        if type(phases) is not list:  # Ensure it's a list
            phases = [phases]

        # Initialize pore and throat data dictionary with conns and coords
        pdata = {}
        tdata = {}

        # Gather list of prop names from network and geometries
        pprops = set(network.props(element='pore', mode=['all', 'deep']) +
                     network.labels(element='pore'))
        tprops = set(network.props(element='throat', mode=['all', 'deep']) +
                     network.labels(element='throat'))

        # Select data from network and geometries using keys
        for item in pprops:
            pdata.update({item: network[item]})
        for item in tprops:
            tdata.update({item: network[item]})

        # Gather list of prop names from phases and physics
        for phase in phases:
            # Gather list of prop names
            pprops = set(phase.props(element='pore', mode=['all', 'deep']) +
                         phase.labels(element='pore'))
            tprops = set(phase.props(element='throat', mode=['all', 'deep']) +
                         phase.labels(element='throat'))
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


class CSV(GenericIO):
    r"""
    This class is used for reading and writing CSV files containing pore and
    throat property data.  This class uses Pandas for transferring data from
    the OpenPNM format to CSV.

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

    @classmethod
    def save(cls, network, filename='', phases=[]):
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
        if type(phases) is not list:  # Ensure it's a list
            phases = [phases]

        dataframes = Pandas.get_data_frames(network=network, phases=phases)
        dfp = dataframes['pore.DataFrame']
        dft = dataframes['throat.DataFrame']
        b = dft.join(other=dfp, how='left')

        # Write to file
        if filename == '':
            filename = network.name
        with cls._write_file(filename=filename, ext='csv') as f:
            b.to_csv(f, index=False)

    @classmethod
    def load(cls, filename, network=None, return_geometry=False):
        r"""
        Opens a 'csv' file, reads in the data, and adds it to the **Network**

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

        with cls._read_file(filename=filename, ext='csv') as f:
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

        if network is None:
            network = OpenPNM.Network.GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)
        return network


class NetworkX(GenericIO):
    r"""
    This class is meant specifcally for exchanging data with NetworkX, which
    is a common tool for dealing with network structures.  A network object
    in NetworkX has a ``to_yaml`` method which produces the correct file format
    for use here.

    Notes
    -----
    1. Each node in a NetworkX object (i.e. ``net``) can be assigned properties
    using syntax like ``net.node[n]['diameter'] = 0.5`` where ``n`` is the
    node number.  There is no need to precede the property name with any
    indication that it is pore data such as \'pore\_\'.  OpenPNM will prepend
    \'pore.\' to each property name.

    2. Since \'pore.coords\' is so central to OpenPNM it should be specified
    in the NetworkX object as \'coords\', and the [X, Y, Z] coordinates of
    each node should be a 3x1 list.

    3. Edges in a NetworkX object are accessed using the index numbers of the
    two nodes it connects, such as ``net.edge[2][3]['length'] = 0.1``
    indicating the edge that connects nodes 2 and 3.  There is no need to
    precede the property name with any indication that it is throat data such
    as \'throat\_\'.  OpenPNM will prepend \'throat.\' to each property name.

    4. The \'throat.conns\' property is essential to OpenPNM, but this does NOT
    need to be specified explicitly as a property in NetworkX.  The
    connectivity is embedded into the network representation in the 'yaml' file
    and is extracted by OpenPNM.
    """

    @classmethod
    def load(cls, filename, network=None, return_geometry=False):
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

        # Open file and read first line, to prevent NetworkX instantiation
        with cls._read_file(filename=filename, ext='yaml') as f:
            line = f.readline()
            if line.startswith('!!python/object:networkx.classes.graph.Graph'):
                a = _yaml.safe_load(f)
            else:
                raise ('Provided file does not appear to be a NetworkX file')

        # Parsing node data
        Np = len(a['node'])
        net.update({'pore.all': _sp.ones((Np,), dtype=bool)})
        for n in a['node'].keys():
            props = a['node'][n]
            for item in props.keys():
                val = a['node'][n][item]
                dtype = type(val)
                # Remove prepended pore. and pore_ if present
                for b in ['pore.', 'pore_']:
                    item = item.replace(b, '')
                # Create arrays for subsequent indexing, if not present already
                if 'pore.'+item not in net.keys():
                    if dtype is list:
                        dtype = type(val[0])
                        cols = len(val)
                        net['pore.'+item] = _sp.ndarray((Np, cols), dtype=dtype)
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
                dtype = type(val)
                # Remove prepended throat. and throat_ if present
                for b in ['throat.', 'throat_']:
                    item = item.replace(b, '')
                # Create arrays for subsequent indexing, if not present already
                if 'throat.'+item not in net.keys():
                    if dtype is list:
                        dtype = type(val[0])
                        cols = len(val)
                        net['throat.'+item] = _sp.ndarray((Nt, cols), dtype=dtype)
                    else:
                        net['throat.'+item] = _sp.ndarray((Nt,), dtype=dtype)
                net['throat.'+item][i] = val
            i += 1

        if network is None:
            network = OpenPNM.Network.GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)
        return network


class iMorph(GenericIO):
    r"""
    Combines two output files from the iMorph program to build a pore network.
    throats_cellsThroatsGraph_Nodes.txt - stores node shape and type information
    throats_cellsThroatsGraph.txt - stores node connectivity
    """

    @classmethod
    def load(cls, path,
             node_file="throats_cellsThroatsGraph_Nodes.txt",
             graph_file="throats_cellsThroatsGraph.txt",
             network=None, voxel_size=None, return_geometry=False):
        r"""
        Loads network data from an iMorph processed image stack

        Parameters
        ----------
        path : string
            The path of the folder where the subfiles are held

        node_file : string
            The file that describes the pores and throats, the
            default iMorph name is: throats_cellsThroatsGraph_Nodes.txt

        graph_file : string
            The file that describes the connectivity of the network, the
            default iMorph name is: throats_cellsThroatsGraph.txt

        network : OpenPNM Network Object
            The OpenPNM Network onto which the data should be loaded.  If no
            network is supplied then an empty import network is created and
            returned.

        voxel_size : float
            Allows the user to define a voxel size different than what is
            contained in the node_file. The value must be in meters.

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
        #
        node_file = _os.path.join(path, node_file)
        graph_file = _os.path.join(path, graph_file)
        # parsing the nodes file
        with open(node_file, 'r') as file:
            Np = _sp.fromstring(file.readline().rsplit('=')[1], sep='\t',
                                dtype=int)[0]
            vox_size = _sp.fromstring(file.readline().rsplit(')')[1], sep='\t',)[0]
            #
            # network always recreated to prevent errors
            network = OpenPNM.Network.Empty(Np=Np, Nt=0)
            #
            # Define expected properies
            network['pore.volume'] = _sp.nan
            scrap_lines = [file.readline() for line in range(4)]
            while True:
                vals = file.readline().split('\t')
                if len(vals) == 1:
                    break
                network['pore.volume'][int(vals[0])] = float(vals[3])
                if 'pore.'+vals[2] not in network.labels():
                    network['pore.'+vals[2]] = False
                network['pore.'+vals[2]][int(vals[0])] = True

        if voxel_size is None:
            voxel_size = vox_size * 1.0E-6  # file stores value in microns

        if voxel_size < 0:
            raise(Exception('Error - Voxel size must be specfied in ' +
                            'the Nodes file or as a keyword argument.'))

        # parsing the graph file
        with open(graph_file, 'r') as file:
            # Define expected properties
            network['pore.coords'] = _sp.zeros((Np, 3))*_sp.nan
            network['pore.types'] = _sp.nan
            network['pore.color'] = _sp.nan
            network['pore.radius'] = _sp.nan
            network['pore.dmax'] = _sp.nan
            network['pore.node_number'] = _sp.nan
            # Scan file to get pore coordinate data
            scrap_lines = [file.readline() for line in range(3)]
            line = file.readline()
            xmax = 0.0
            ymax = 0.0
            zmax = 0.0
            node_num = 0
            while line != 'connectivity table\n':
                vals = _sp.fromstring(line, sep='\t')
                xmax = vals[1] if vals[1] > xmax else xmax
                ymax = vals[2] if vals[2] > ymax else ymax
                zmax = vals[3] if vals[3] > zmax else zmax
                network['pore.coords'][int(vals[0]), :] = vals[1:4]
                network['pore.types'][int(vals[0])] = vals[4]
                network['pore.color'][int(vals[0])] = vals[5]
                network['pore.radius'][int(vals[0])] = vals[6]
                network['pore.dmax'][int(vals[0])] = vals[7]
                network['pore.node_number'][int(vals[0])] = node_num
                node_num += 1
                line = file.readline()
            # Scan file to get to connectivity data
            scrap_lines.append(file.readline())  # Skip line
            # Create sparse lil array incrementally build adjacency matrix
            lil = _sp.sparse.lil_matrix((Np, Np), dtype=int)
            while True:
                vals = _sp.fromstring(file.readline(), sep='\t', dtype=int)
                if len(vals) <= 1:
                    break
                lil.rows[vals[0]] = vals[2:]
                lil.data[vals[0]] = _sp.ones(vals[1])

        # fixing any negative volumes or distances so they are 1 voxel/micron
        network['pore.volume'][_sp.where(network['pore.volume'] < 0)[0]] = 1.0
        network['pore.radius'][_sp.where(network['pore.radius'] < 0)[0]] = 1.0
        network['pore.dmax'][_sp.where(network['pore.dmax'] < 0)[0]] = 1.0

        # Add adjacency matrix to OpenPNM network
        conns = _sp.sparse.triu(lil, k=1, format='coo')
        network.update({'throat.all': _sp.ones(len(conns.col), dtype=bool)})
        network['throat.conns'] = _sp.vstack([conns.row, conns.col]).T

        network['pore.to_trim'] = False
        network['pore.to_trim'][network.pores('*throat')] = True
        Ts = network.pores('to_trim')
        new_conns = network.find_neighbor_pores(pores=Ts, flatten=False)
        network.extend(throat_conns=new_conns, labels='new_conns')
        for item in network.props('pore'):
            item = item.split('.')[1]
            arr = _sp.ones_like(network['pore.'+item])[0]
            arr = _sp.tile(A=arr, reps=[network.Nt, 1])*_sp.nan
            network['throat.'+item] = _sp.squeeze(arr)
            network['throat.'+item][network.throats('new_conns')] = \
                network['pore.'+item][Ts]
        network.trim(pores=Ts)

        # setting up boundary pores
        x_coord, y_coord, z_coord = _sp.hsplit(network['pore.coords'], 3)
        network['pore.front_boundary'] = _sp.ravel(x_coord == 0)
        network['pore.back_boundary'] = _sp.ravel(x_coord == xmax)
        network['pore.left_boundary'] = _sp.ravel(y_coord == 0)
        network['pore.right_boundary'] = _sp.ravel(y_coord == ymax)
        network['pore.bottom_boundary'] = _sp.ravel(z_coord == 0)
        network['pore.top_boundary'] = _sp.ravel(z_coord == zmax)

        # removing any pores that got classified as a boundary pore but
        # weren't labled a border_cell_face
        ps = _sp.where(~_sp.in1d(network.pores('*_boundary'),
                                 network.pores('border_cell_face')))[0]
        ps = network.pores('*_boundary')[ps]
        for side in ['front', 'back', 'left', 'right', 'top', 'bottom']:
            network['pore.'+side+'_boundary'][ps] = False
        # setting internal label
        network['pore.internal'] = False
        network['pore.internal'][network.pores('*_boundary', mode='not')] = True

        # adding props to border cell face throats and from pores
        Ts = _sp.where(network['throat.conns'][:, 1] >
                       network.pores('border_cell_face')[0] - 1)[0]
        faces = network['throat.conns'][Ts, 1]
        for item in network.props('pore'):
            item = item.split('.')[1]
            network['throat.'+item][Ts] = network['pore.'+item][faces]
        network['pore.volume'][faces] = 0.0

        # applying unit conversions
        # TODO: Determine if radius and dmax are indeed microns and not voxels
        network['pore.coords'] = network['pore.coords'] * 1e-6
        network['pore.radius'] = network['pore.radius'] * 1e-6
        network['pore.dmax'] = network['pore.dmax'] * 1e-6
        network['pore.volume'] = network['pore.volume'] * voxel_size**3
        network['throat.coords'] = network['throat.coords'] * 1e-6
        network['throat.radius'] = network['throat.radius'] * 1e-6
        network['throat.dmax'] = network['throat.dmax'] * 1e-6
        network['throat.volume'] = network['throat.volume'] * voxel_size**3
        #
        # checking network health to generate warnings for the user
        network.health_dict = network.check_network_health()
        logger.info('Network health stored as network.health_dict')
        if return_geometry:
            geometry = cls.fetch_geometry(network)
            network = (network, geometry)
        return network


class MARock(GenericIO):
    r"""
    3DMA-Rock is a network extraction algorithm developed by Brent Lindquist
    and his group [1].  It uses Medial Axis thinning to find the skeleton of
    the pore space, then extracts geometrical features such as pore volume and
    throat cross-sectional area.

    [1] Lindquist, W. Brent, S. M. Lee, W. Oh, A. B. Venkatarangan, H. Shin,
    and M. Prodanovic. "3DMA-Rock: A software package for automated analysis
    of rock pore structure in 3-D computed microtomography images." SUNY Stony
    Brook (2005).
    """

    @classmethod
    def load(cls, path, network=None, voxel_size=1, return_geometry=False):
        r"""
        Load data from a 3DMA-Rock extracted network.  This format consists of
        two files: 'rockname.np2th' and 'rockname.th2pn'.  They should be
        stored together in a folder which is referred to by the path argument.
        These files are binary and therefore not human readable.

        Parameters
        ----------
        path : string
            The location of the 'np2th' and 'th2np' files. This can be an
            absolute path or relative to the current working directory.

        network : OpenPNM Network Object
            If an Network object is recieved, this method will add new data to
            it but NOT overwrite anything that already exists.  This can be
            used to append data from different sources.

        voxel_size : scalar
            The resolution of the image on which 3DMA-Rock was run, in terms of
            the linear length of eac voxel. The default is 1.  This is used to
            scale the voxel counts to actual dimension. It is recommended that
            this value be in SI units [m] to work well with OpenPNM.

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

        for file in _os.listdir(path):
            if file.endswith(".np2th"):
                np2th_file = _os.path.join(path, file)
            elif file.endswith(".th2np"):
                th2np_file = _os.path.join(path, file)

        with open(np2th_file, mode='rb') as f:
            [Np, Nt] = _sp.fromfile(file=f, count=2, dtype='u4')
            net['pore.boundary_type'] = _sp.ndarray([Np, ], int)
            net['throat.conns'] = _sp.ones([Nt, 2], int)*(-1)
            net['pore.coordination'] = _sp.ndarray([Np, ], int)
            net['pore.ID_number'] = _sp.ndarray([Np, ], int)
            for i in range(0, Np):
                ID = _sp.fromfile(file=f, count=1, dtype='u4')
                net['pore.ID_number'][i] = ID
                net['pore.boundary_type'][i] = _sp.fromfile(file=f, count=1,
                                                            dtype='u1')
                z = _sp.fromfile(file=f, count=1, dtype='u4')[0]
                net['pore.coordination'][i] = z
                att_pores = _sp.fromfile(file=f, count=z, dtype='u4')
                att_throats = _sp.fromfile(file=f, count=z, dtype='u4')
                for j in range(0, len(att_throats)):
                    t = att_throats[j] - 1
                    p = att_pores[j] - 1
                    net['throat.conns'][t] = [i, p]
            net['throat.conns'] = _sp.sort(net['throat.conns'], axis=1)
            net['pore.volume'] = _sp.fromfile(file=f, count=Np, dtype='u4')
            nx = _sp.fromfile(file=f, count=1, dtype='u4')
            nxy = _sp.fromfile(file=f, count=1, dtype='u4')
            pos = _sp.fromfile(file=f, count=Np, dtype='u4')
            ny = nxy/nx
            ni = _sp.mod(pos, nx)
            nj = _sp.mod(_sp.floor(pos/nx), ny)
            nk = _sp.floor(_sp.floor(pos/nx)/ny)
            net['pore.coords'] = _sp.array([ni, nj, nk]).T

        with open(th2np_file, mode='rb') as f:
            Nt = _sp.fromfile(file=f, count=1, dtype='u4')[0]
            net['throat.area'] = _sp.ones([Nt, ], dtype=int)*(-1)
            for i in range(0, Nt):
                ID = _sp.fromfile(file=f, count=1, dtype='u4')
                net['throat.area'][i] = _sp.fromfile(file=f, count=1,
                                                     dtype='f4')
                numvox = _sp.fromfile(file=f, count=1, dtype='u4')
                att_pores = _sp.fromfile(file=f, count=2, dtype='u4')
            nx = _sp.fromfile(file=f, count=1, dtype='u4')
            nxy = _sp.fromfile(file=f, count=1, dtype='u4')
            pos = _sp.fromfile(file=f, count=Nt, dtype='u4')
            ny = nxy/nx
            ni = _sp.mod(pos, nx)
            nj = _sp.mod(_sp.floor(pos/nx), ny)
            nk = _sp.floor(_sp.floor(pos/nx)/ny)
            net['throat.coords'] = _sp.array([ni, nj, nk]).T
            net['pore.internal'] = net['pore.boundary_type'] == 0

        # Convert voxel area and volume to actual dimensions
        net['throat.area'] = (voxel_size**2)*net['throat.area']
        net['pore.volume'] = (voxel_size**3)*net['pore.volume']

        if network is None:
            network = OpenPNM.Network.GenericNetwork()
        network = cls._update_network(network=network, net=net,
                                      return_geometry=return_geometry)

        # Trim headless throats before returning
        ind = _sp.where(network['throat.conns'][:, 0] == -1)[0]
        network.trim(throats=ind)

        return network
