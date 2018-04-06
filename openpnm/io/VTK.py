from xml.etree import ElementTree as ET
import numpy as np
from openpnm.core import logging, Workspace
from openpnm.network import GenericNetwork
from openpnm.io import GenericIO, Dict
from openpnm.utils import FlatDict
logger = logging.getLogger(__name__)
ws = Workspace()


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
    def save(cls, network, phases=[], filename='', delim=' | '):
        r"""
        Save network and phase data to a single vtp file for visualizing in
        Paraview

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network containing the data to be written

        phases : list, optional
            A list containing OpenPNM Phase object(s) containing data to be
            written

        filename : string, optional
            Filename to write data.  If no name is given the file is named
            after ther network

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)

        am = Dict.to_dict(network=network, phases=phases, interleave=True,
                          categorize_by=['object', 'data'])
        am = FlatDict(am, delimiter=delim)
        key_list = list(sorted(am.keys()))

        network = network[0]
        points = network['pore.coords']
        pairs = network['throat.conns']
        num_points = np.shape(points)[0]
        num_throats = np.shape(pairs)[0]

        root = ET.fromstring(VTK._TEMPLATE)
        piece_node = root.find('PolyData').find('Piece')
        piece_node.set("NumberOfPoints", str(num_points))
        piece_node.set("NumberOfLines", str(num_throats))
        points_node = piece_node.find('Points')
        coords = VTK._array_to_element("coords", points.T.ravel('F'), n=3)
        points_node.append(coords)
        lines_node = piece_node.find('Lines')
        connectivity = VTK._array_to_element("connectivity", pairs)
        lines_node.append(connectivity)
        offsets = VTK._array_to_element("offsets", 2*np.arange(len(pairs))+2)
        lines_node.append(offsets)

        point_data_node = piece_node.find('PointData')
        for key in key_list:
            array = am[key]
            if array.dtype == np.bool:
                array = array.astype(int)
            if array.size != num_points:
                continue
            element = VTK._array_to_element(key, array)
            if element is not None:
                point_data_node.append(element)

        cell_data_node = piece_node.find('CellData')
        for key in key_list:
            array = am[key]
            if array.dtype == np.bool:
                array = array.astype(int)
            if array.size != num_throats:
                continue
            element = VTK._array_to_element(key, array)
            if element is not None:
                cell_data_node.append(element)

        if filename == '':
            filename = project.name
        filename = cls._parse_filename(filename=filename, ext='vtp')

        tree = ET.ElementTree(root)
        tree.write(filename)

        with open(filename, 'r+') as f:
            string = f.read()
            string = string.replace('</DataArray>', '</DataArray>\n\t\t\t')
            f.seek(0)
            # consider adding header: '<?xml version="1.0"?>\n'+
            f.write(string)

    @classmethod
    def load(cls, filename, project=None, delim=' | '):
        r"""
        Read in pore and throat data from a saved VTK file.

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project is supplied then one will be created and returned.

        """
        net = {}

        filename = cls._parse_filename(filename, ext='vtp')
        tree = ET.parse(filename)
        piece_node = tree.find('PolyData').find('Piece')

        # Extract connectivity
        conn_element = piece_node.find('Lines').find('DataArray')
        conns = VTK._element_to_array(conn_element, 2)
        # Extract coordinates
        coord_element = piece_node.find('Points').find('DataArray')
        coords = VTK._element_to_array(coord_element, 3)

        # Extract pore data
        for item in piece_node.find('PointData').iter('DataArray'):
            key = item.get('Name')
            array = VTK._element_to_array(item)
            net[key] = array
        # Extract throat data
        for item in piece_node.find('CellData').iter('DataArray'):
            key = item.get('Name')
            array = VTK._element_to_array(item)
            net[key] = array

        if project is None:
            project = ws.new_project()
        project = Dict.from_dict(dct=net, project=project, delim=delim)

        # Clean up data values, if necessary, like convert array's of
        # 1's and 0's into boolean.
        project = cls._convert_data(project)

        # Add coords and conns to network
        network = project.network
        network.update({'throat.conns': conns})
        network.update({'pore.coords': coords})

        return project

    @classmethod
    def _array_to_element(cls, name, array, n=1):
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
        element = None
        if str(array.dtype) in dtype_map.keys():
            element = ET.Element('DataArray')
            element.set("Name", name)
            element.set("NumberOfComponents", str(n))
            element.set("type", dtype_map[str(array.dtype)])
            element.text = '\t'.join(map(str, array.ravel()))
        return element

    @classmethod
    def _element_to_array(cls, element, n=1):
        string = element.text
        dtype = element.get("type")
        array = np.fromstring(string, sep='\t')
        array = array.astype(dtype)
        if n is not 1:
            array = array.reshape(array.size//n, n)
        return array
