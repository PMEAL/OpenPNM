from xml.etree import ElementTree as _ET
import OpenPNM.Utilities.misc as misc
import numpy as _np

class VTK():
    r"""
    Class for writing a Vtp file to be read by ParaView

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network containing the data to be written

    filename : string, optional
        Filename to write data.  If no name is given the file is named after
        ther network

    phase : list, optional
        A list contain OpenPNM Phase object(s) containing data to be written

    """

    def __init__(self,network,filename='',phases=[],**kwargs):
        r"""
        Initialize
        """
        self._TEMPLATE = '''
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
        if filename == '':
            filename = network.name+'.vtp'
        self._net = network
        self._phases = phases
        self._write(filename)

    def _array_to_element(self, name, array, n=1):
        dtype_map = {
            'int8'   : 'Int8',
            'int16'  : 'Int16',
            'int32'  : 'Int32',
            'int64'  : 'Int64',
            'uint8'  : 'UInt8',
            'uint16' : 'UInt16',
            'uint32' : 'UInt32',
            'uint64' : 'UInt64',
            'float32': 'Float32',
            'float64': 'Float64',
            'str'    : 'String',
        }
        element = _ET.Element('DataArray')
        element.set("Name", name)
        element.set("NumberOfComponents", str(n))
        element.set("type", dtype_map[str(array.dtype)])
        element.text = '\t'.join(map(str,array.ravel()))
        return element

    def _element_to_array(self, element, n=1):
        string = element.text
        dtype = element.get("type")
        array = _np.fromstring(string, sep='\t')
        array = array.astype(dtype)
        if n is not 1:
            array = array.reshape(array.size//n, n)
        return array

    def _write(self,filename):
        r"""
        Write Network to a VTK file for visualizing in Paraview

        Parameters
        ----------

        network : OpenPNM Network Object

        filename : string
            Full path to desired file location

        phases : Phases that have properties we want to write to file

        """
        phases = self._phases
        network = self._net

        root = _ET.fromstring(self._TEMPLATE)
        objs = []
        if _np.shape(phases)==():
            phases = [phases]
        for phase in phases:
            objs.append(phase)
        objs.append(network)
        am = misc.amalgamate_data(objs=objs)
        key_list = list(sorted(am.keys()))
        points = am[network.name+'.pore.coords']
        pairs = network['throat.conns']

        num_points = len(points)
        num_throats = len(pairs)

        piece_node = root.find('PolyData').find('Piece')
        piece_node.set("NumberOfPoints", str(num_points))
        piece_node.set("NumberOfLines", str(num_throats))

        points_node = piece_node.find('Points')
        coords = self._array_to_element("coords", points.T.ravel('F'), n=3)
        points_node.append(coords)

        lines_node = piece_node.find('Lines')
        connectivity = self._array_to_element("connectivity", pairs)
        lines_node.append(connectivity)
        offsets = self._array_to_element("offsets", 2*_np.arange(len(pairs))+2)
        lines_node.append(offsets)

        point_data_node = piece_node.find('PointData')
        for key in key_list:
            array = am[key]
            if array.dtype == _np.bool: array = array.astype(int)
            if array.size != num_points: continue
            element = self._array_to_element(key, array)
            point_data_node.append(element)

        cell_data_node = piece_node.find('CellData')
        for key in key_list:
            array = am[key]
            if array.dtype == _np.bool: array = array.astype(int)
            if array.size != num_throats: continue
            element = self._array_to_element(key, array)
            cell_data_node.append(element)

        tree = _ET.ElementTree(root)
        tree.write(filename)

        #Make pretty
        with open(filename, "r+") as f:
            string = f.read()
            string = string.replace("</DataArray>", "</DataArray>\n\t\t\t")
            f.seek(0)
            # consider adding header: '<?xml version="1.0"?>\n'+
            f.write(string)

    def read(self,filename):
        r'''
        Read in pore and throat data from a saved VTK file.

        Notes
        -----
        This will NOT reproduce original simulation, since all models and object
        relationships are lost.  Use IO.Save and IO.Load for that.
        '''
        network = {}
        tree = _ET.parse(filename)
        piece_node = tree.find('PolyData').find('Piece')

        # extract connectivity
        conn_element = piece_node.find('Lines').find('DataArray')
        array = self._element_to_array(conn_element, 2)
        network['heads'], network['tails'] = array.T

        for element in piece_node.find('PointData').iter('DataArray'):

            key = element.get('Name')
            array = self._element_to_array(element)
            network[key] = array

        return network