from xml.etree import ElementTree as ET

import numpy as np
from OpenPNM.Visualization import GenericVisualization


class VTK(GenericVisualization):
    r"""
    VTK - Class for writing a Vtp file to be read by ParaView

    
    Examples
    --------
    Create and store a basic network.

    >>> import OpenPNM as PNM
    >>> net = PNM.Generators.SimpleCubic(divisions = [40,40,20],shape=[0.,4.,0.,4.,0.,2.]).generate()
    >>> vis = PNM.Visualization.VTK()
    >>> vis.write(net)

    """

    def __init__(self,**kwargs):
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
                </Piece>
            </PolyData>
        </VTKFile>
        '''.strip()
        
        super(VTK,self).__init__(**kwargs)
#        self._logger.debug("Execute constructor")
    
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
        element = ET.Element('DataArray')
        element.set("Name", name)
        element.set("NumberOfComponents", str(n))
        element.set("type", dtype_map[str(array.dtype)])
        element.text = '\t'.join(map(str,array.ravel()))
        return element
    
    def _element_to_array(self, element, n=1):
        string = element.text
        dtype = element.get("type")
        array = np.fromstring(string, sep='\t')
        array = array.astype(dtype)
        if n is not 1:
            array = array.reshape(array.size//n, n)
        return array
    
    def write(self, network, filename='output_file.vtp', fluids=[], pretty=True):
        r"""
        Write Network to a VTK file for visualizing in Paraview
    
        Parameters
        ----------
    
        network : OpenPNM Network Object
    
        filename : string
            Full path to desired file location
            
        fluids : Fluids that have properties we want to write to file
    
        pretty : Add linebreaks at the end of tag closures
        """
        
        root = ET.fromstring(self._TEMPLATE)
        objs = []
        if np.shape(fluids)==():
            fluids = [fluids]
        for fluid in fluids:
            objs.append(fluid)
        objs.append(network)
        am = network.amalgamate_data(objs=objs)
        key_list = list(sorted(am.keys()))
        points = am[network.name+'.pore.coords']
        pairs = network.get_throat_data(prop='conns')
    
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
        offsets = self._array_to_element("offsets", 2*np.arange(len(pairs))+2)
        lines_node.append(offsets)
    
        point_data_node = piece_node.find('PointData')
        for key in key_list:
            array = am[key]
            if array.dtype == np.bool: array = array.astype(int)
            if array.size != num_points: continue
            element = self._array_to_element(key, array)
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
        array = self._element_to_array(conn_element, 2)
        network['heads'], network['tails'] = array.T
    
        for element in piece_node.find('PointData').iter('DataArray'):
    
            key = element.get('Name')
            array = self._element_to_array(element)
            network[key] = array
    
        return network