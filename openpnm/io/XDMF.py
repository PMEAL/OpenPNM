import h5py
import xml.etree.cElementTree as ET
from flatdict import FlatDict
from openpnm.io import Dict, GenericIO


class XDMF(GenericIO):
    r"""
    The eXtensible Data Model Format combines XML descriptors with HDF5 data
    storage


    For more information visit the webiste:
    `XDMF.org <http://www.xdmf.org/index.php/Main_Page>`_
    """

    _header = '''<?xml version="1.0" ?>
                 <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'''

    @classmethod
    def save(cls, network, phases=[], filename=''):
        r"""
        Saves data from the given objects into the specified file.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)
            A list of phase objects whose data are to be included

        Notes
        -----
        This method only saves the data, not any of the pore-scale models or
        other attributes.  To save an actual OpenPNM Project use the
        ``Workspace`` object.

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        network = network[0]

        if filename == '':
            filename = project.name
        path = cls._parse_filename(filename=filename, ext='xdf')
        fname_xdf = path.name
        fname_hdf = path.stem+".hdf"
        path = path.parent
        f = h5py.File(path.joinpath(fname_hdf), "w")

        d = Dict.to_dict(network, phases=phases, interleave=True,
                         flatten=False, categorize_by=['element', 'data'])

        # Make HDF5 file with all datasets, and no groups
        D = FlatDict(d, delimiter='/')
        for item in D.keys():
            if 'U' in str(D[item][0].dtype):
                pass
            else:
                f.create_dataset(name='/'+item, shape=D[item].shape,
                                 dtype=D[item].dtype, data=D[item])
        # Add coordinate and connection information to top of HDF5 file
        f["coordinates"] = network["pore.coords"]
        f["connections"] = network["throat.conns"]

        # setup xdmf file
        root = create_root('Xdmf')
        domain = create_domain()
        grid = create_grid(Name="Structure", GridType="Uniform")

        # geometry coordinates
        row, col = f["coordinates"].shape
        dims = ' '.join((str(row), str(col)))
        hdf_loc = f.filename + ":coordinates"
        geo_data = create_data_item(value=hdf_loc, Dimensions=dims,
                                    Format='HDF', DataType="Float")
        geo = create_geometry(GeometryType="XYZ")
        geo.append(geo_data)

        # topolgy connections
        row, col = f["connections"].shape  # col first then row
        dims = ' '.join((str(row), str(col)))
        hdf_loc = f.filename + ":connections"
        topo_data = create_data_item(value=hdf_loc, Dimensions=dims,
                                     Format="HDF", NumberType="Int")
        topo = create_topology(TopologyType="Polyline",
                               NodesPerElement=str(2),
                               NumberOfElements=str(row))
        topo.append(topo_data)

        # Add pore and throat properties
        for item in D.keys():
            if item not in ['coordinates', 'connections']:
                attr_type = 'Scalar'
                shape = f[item].shape
                dims = ''.join([str(i) + ' ' for i in list(shape)[::-1]])
                hdf_loc = fname_hdf + ":" + item
                attr = create_data_item(value=hdf_loc,
                                        Dimensions=dims,
                                        Format='HDF',
                                        Precision='8',
                                        DataType='Float')
                name = item.replace('/', ' | ')
                if 'throat' in item:
                    Center = "Cell"
                else:
                    Center = "Node"
                el_attr = create_attribute(Name=name, Center=Center,
                                           AttributeType=attr_type)
                el_attr.append(attr)
                grid.append(el_attr)

        grid.append(topo)
        grid.append(geo)
        domain.append(grid)
        root.append(domain)

        with open(path.joinpath(fname_xdf), 'w') as file:
            file.write(cls._header)
            file.write(ET.tostring(root).decode("utf-8"))

        # CLose the HDF5 file
        f.close()


def create_root(Name):
    return ET.Element(Name)


def create_domain():
    return ET.Element("Domain")


def create_geometry(GeometryType, **attribs):
    element = ET.Element('Geometry')
    element.attrib.update({'GeometryType': GeometryType})
    element.attrib.update(attribs)
    return element


def create_topology(TopologyType, NumberOfElements, **attribs):
    element = ET.Element('Topology')
    element.attrib.update({'TopologyType': TopologyType,
                           'NumberOfElements': NumberOfElements})
    if TopologyType in ['Polyline']:
        if 'NodesPerElement' not in attribs.keys():
            raise Exception('NodesPerElement must be specified')
    element.attrib.update(attribs)
    return element


def create_attribute(Name, **attribs):
    element = ET.Element('Attribute')
    element.attrib.update({'Name': Name})
    element.attrib.update(attribs)
    return element


def create_time(type='Single', Value=None):
    element = ET.Element('Attribute')
    if type == 'Single' and Value:
        element.attrib['Value'] = Value
    return element


def create_grid(Name, GridType, **attribs):
    element = ET.Element('Grid')
    element.attrib.update({'Name': Name,
                           'GridType': GridType,
                           'Section': None})
    element.attrib.update(attribs)
    if element.attrib['GridType'] is not 'Subset':
        if 'Section' in element.attrib.keys():
            del element.attrib['Section']
    return element


def create_data_item(value, Dimensions, **attribs):
    element = ET.Element('DataItem')
    element.attrib.update({'ItemType': "Uniform",
                           'Format': "XML",
                           'DatarType': "Float",
                           'Precision': "4",
                           'Rank': "1",
                           'Dimensions': Dimensions,
                           'Function': None})
    element.attrib.update(attribs)
    if element.attrib['Function'] is None:
        del element.attrib['Function']
    element.text = value
    return element
