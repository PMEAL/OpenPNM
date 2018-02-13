import xml.etree.cElementTree as ET
from openpnm.io import HDF5


class XDMF:

    header = '''<?xml version="1.0" ?>
                <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'''

    def __init__(self, Version="2.0"):
        self.Version = Version

    @classmethod
    def save(cls, network, phases=[], filename=''):
        f = HDF5.to_hdf5(network, phases=phases, filename=filename,
                         interleave=True, flatten=False,
                         categorize_objects=False, categorize_data=False)
        filename = f.filename.split('.')[0]

        # Add coordinate and connection information to top of HDF5 file
        f["coordinates/pore"] = network["pore.coords"]
        f["connections/throat"] = network["throat.conns"]

        # setup xdmf file
        root = ET.Element("xdmf")
        domain = create_domain()
        grid = create_grid(Name="Structure", GridType="Uniform")

        # geometry coordinates
        col, row = f["coordinates/pore"].shape
        dims = create_dimension(network["pore.coords"].shape)
        hdf_loc = f.filename + ":coordinates/pore"
        geo_data = create_data_item(value=hdf_loc, Dimensions=dims,
                                    Format='HDF', NumberType="Float")
        geo = create_geometry(GeometryType="XYZ")
        geo.append(geo_data)

        # topolgy connections
        row, col = f["connections/throat"].shape  # col first then row
        dims = str(row) + " " + str(col)
        hdf_loc = f.filename + ":connections/throat"
        top_data = create_data_item(value=hdf_loc, Dimensions=dims,
                                    Format="HDF", NumberType="Int")
        topo = create_topology(TopologyType="Polyline",
                               NumberOfElements=str(row))
        topo.append(top_data)

        # attributes
#        for key in network.labels():
#            shape = f["labels/" + "/".join(key.split("."))].shape
#            dimensions = create_dimension(shape)
#            hdf_loc = network.name + ".h5:labels/" + "/".join(key.split("."))
#            scalar = create_data_item(value=hdf_loc,
#                                      Dimensions=str(dimensions),
#                                      Format='HDF')
#            el_attr = create_attribute(Name=key)
#            el_attr.append(scalar)
#            grid.append(el_attr)

        # properties, same method as create Prop
        for obj in f.keys():
            for propname in f[obj].keys():
                shape = f[obj + '/' + propname].shape
                dimensions = create_dimension(shape)
                hdf_loc = f.filename + ":" + obj + '/' + propname
                scalar = create_data_item(value=hdf_loc,
                                          Dimensions=str(dimensions),
                                          Format='HDF',
                                          Precision='8',
                                          NumberType='Float')
                el_attr = create_attribute(Name=propname)
                el_attr.append(scalar)
                grid.append(el_attr)

        grid.append(topo)
        grid.append(geo)
        domain.append(grid)
        root.append(domain)

        with open(filename+'.xdmf', 'w') as file:
            file.write(cls.header)
            file.write(ET.tostring(root).decode("utf-8"))


def create_dimension(shape):
    return ''.join([str(i) + ' ' for i in list(shape)[::-1]])


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
    element.attrib.update(attribs)
    return element


def create_attribute(Name, **attribs):
    element = ET.Element('Attribute')
    element.attrib.update({'Name': Name,
                           'AttributeType': 'Scalar',
                           'Center': 'Node'})
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
                           'NumberType': "Float",
                           'Precision': "4",
                           'Rank': "1",
                           'Dimensions': Dimensions,
                           'Function': None})
    element.attrib.update(attribs)
    if element.attrib['Function'] is None:
        del element.attrib['Function']
    element.text = value
    return element
