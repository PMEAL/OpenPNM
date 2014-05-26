from xml.etree import ElementTree as ET

import numpy as np

TEMPLATE = '''
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
        </Piece>
    </PolyData>
</VTKFile>
'''.strip()

def _array_to_element(name, array, n=1):
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
    element = ET.Element('DataArray')
    element.set("Name", name)
    element.set("NumberOfComponents", str(n))
    element.set("type", dtype_map[str(array.dtype)])
    element.text = '\t'.join(map(str,array.ravel()))
    return element

def _element_to_array(element, n=1):
    string = element.text
    dtype = element.get("type")
    array = np.fromstring(string, sep='\t')
    array = array.astype(dtype)
    if n is not 1:
        array = array.reshape(array.size//n, n)
    return array

def write(filename, network, fluids=[], pretty=True):
    r"""
    Write Network to a VTK file for visualizing in Paraview

    Parameters
    ----------
    filename : string
        Full path to desired file location

    network : OpenPNM Network Object

    fluids : Fluids that have properties we want to write to file

    pretty : Add linebreaks at the end of tag closures
    """
    root = ET.fromstring(TEMPLATE)

    am = network.amalgamate_pore_data(fluids=fluids)
    key_list = list(sorted(am.keys()))
    points = am['pore_coords']
    pairs = network.get_throat_data(prop='conns')

    num_points = len(points)
    num_throats = len(pairs)
    
    piece_node = root.find('PolyData').find('Piece')
    piece_node.set("NumberOfPoints", str(num_points))
    piece_node.set("NumberOfLines", str(num_throats))

    points_node = piece_node.find('Points')
    coords = _array_to_element("coords", points.T.ravel('F'), n=3)
    points_node.append(coords)

    lines_node = piece_node.find('Lines')
    connectivity = _array_to_element("connectivity", pairs)
    lines_node.append(connectivity)
    offsets = _array_to_element("offsets", 2*np.arange(len(pairs))+2)
    lines_node.append(offsets)

    point_data_node = piece_node.find('PointData')
    for key in key_list:
        array = am[key]
        if array.dtype == np.bool: array = array.astype(int)
        if array.size != num_points: continue
        element = _array_to_element(key, array)
        point_data_node.append(element)

    tree = ET.ElementTree(root)
    tree.write(filename)

    if pretty:
        with open(filename, "r+") as f:
            string = f.read()
            string = string.replace("</DataArray>", "</DataArray>\n\t\t\t")
            f.seek(0)
            # consider adding header: '<?xml version="1.0"?>\n'+
            f.write(string)

def read(filename):
    network = {}
    tree = ET.parse(filename)
    piece_node = tree.find('PolyData').find('Piece')

    # extract connectivity
    conn_element = piece_node.find('Lines').find('DataArray')
    array = _element_to_array(conn_element, 2)
    network['heads'], network['tails'] = array.T

    for element in piece_node.find('PointData').iter('DataArray'):

        key = element.get('Name')
        array = _element_to_array(element)
        network[key] = array

    return network