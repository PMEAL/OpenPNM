from OpenPNM.Utilities import misc
import scipy as _sp
import numpy as _np
import os as _os
import pickle as _pickle
from xml.etree import ElementTree as _ET


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
            Filename to write data.  If no name is given the file is named after
            ther network

        phases : list, optional
            A list contain OpenPNM Phase object(s) containing data to be written

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
        >>> geo = OpenPNM.Geometry.Stick_and_Ball(network=pn,
        ...                                       pores=pn.pores(),
        ...                                       throats=pn.throats())
        >>> air = OpenPNM.Phases.Air(network=pn)
        >>> phys = OpenPNM.Physics.Standard(network=pn, phase=air,
        ...                                 pores=pn.pores(), throats=pn.throats())

        >>> import OpenPNM.Utilities.IO as io
        >>> io.VTK.save(pn,'test_pn.vtp',[air])

        >>> # Delete the new file
        >>> import os
        >>> os.remove('test_pn.vtp')
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
        am = misc.amalgamate_data(objs=objs)
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
    def load(filename):
        r"""
        Read in pore and throat data from a saved VTK file.

        Notes
        -----
        This will NOT reproduce original simulation, since all models and object
        relationships are lost.  Use IO.Save and IO.Load for that.
        """
        network = OpenPNM.Network.GenericNetwork()
        tree = _ET.parse(filename)
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
        Write Network to a Mat file for exporting to Matlab. This method will be
        enhanced in a future update, and it's functionality may change!

        Parameters
        ----------

        network : OpenPNM Network Object

        filename : string
            Desired file name, defaults to network name if not given

        phases : list of phase objects ([])
            Phases that have properties we want to write to file

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,
        ...                                     pores=pn.pores(),
        ...                                     throats=pn.throats())
        >>> air = OpenPNM.Phases.TestPhase()
        >>> import OpenPNM.Utilities.IO as io
        >>> io.MAT.save(network=pn, filename='test_pn.mat', phases=air)

        >>> # Remove newly created file
        >>> import os
        >>> os.remove('test_pn.mat')

        """
        if filename == '':
            filename = network.name
        filename = filename.split('.')[0] + '.mat'

        pnMatlab = {}
        new = []
        old = []
        for keys in network.keys():
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

                for keys in phases[j].keys():
                    old.append(keys)
                    new.append(phases[j].name+'_'+keys.replace('.', '_'))

                for i in range(len(phases[j])):
                    pnMatlab[new[i]] = phases[j][old[i]]

        _sp.io.savemat(file_name=filename, mdict=pnMatlab)

    @staticmethod
    def load():
        r"""
        This method is not implemented yet.
        """
        raise NotImplemented()
