import xml.etree.cElementTree as ET
import h5py


class XDMF:
    _START \
        = '''<?xml version="1.0" ?>
    <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
       '''

    def __init__(self, Version="2.0"):
        self.Version = Version

    def get_hdf5(self, network):
        # I would like this to be replaced by a call to the HDF5 class
        f = h5py.File(network.name + ".hdf5", "a")
        for value in network.labels():
            f["labels/" + "/".join(value.split("."))] = 1*network[value]
        # geo, phases, physics
        simulation = network.simulation
        for obj in simulation:
            for value in obj.props():
                arr = value.split(".")
                path = "properties/"+obj.name+"/"+arr[0]+"/"+arr[1]
                f[path] = network[value]
        return f

    def save(self, network, filename):
        f = self.get_hdf5(network=network)

        # Add coordinate and connection information to top of HDF5 file
        f["coordinates/pore"] = network["pore.coords"]
        f["connections/throat"] = network["throat.conns"]

        root = ET.Element("xdmf")
        domain = self.create_domain()
        grid = self.create_grid(Name="Structure", GridType="Uniform")

        # geometry coordinates
        col, row = f["coordinates/pore"].shape
        dims = self.createDimension(network["pore.coords"].shape)
        hdf_loc = network.name + ".h5:coordinates/pore"
        geo_data = self.create_data_item(value=hdf_loc, Dimensions=dims,
                                         Format='HDF')
        geo = self.create_geometry(GeometryType="XYZ")
        geo.append(geo_data)

        # topolgy connections
        row, col = f["connections/throat"].shape  # col first then row
        dims = str(row) + " " + str(col)
        hdf_loc = network.name + ".h5:connections/throat"
        top_data = self.create_data_item(value=hdf_loc, Dimensions=dims,
                                         Format="HDF", NumberType="Int")
        topo = self.create_topology(TopologyType="Polyline",
                                    NumberOfElements=str(row))
        topo.append(top_data)

        # attributes
        for key in network.labels():
            shape = f["labels/" + "/".join(key.split("."))].shape
            dimensions = self.createDimension(shape)
            hdf_loc = network.name + ".h5:labels/" + "/".join(key.split("."))
            scalar = self.create_data_item(value=hdf_loc,
                                           Dimensions=str(dimensions),
                                           Format='HDF')
            el_attr = self.create_attribute(Name=key)
            el_attr.append(scalar)
            grid.append(el_attr)

        # properties, same method as create Prop
        props = [network.geometries, network.phases, network.physics]

        for key in props:
            for k in key:
                obj = key[k]
                for value in obj.models.keys():
                    arr = value.split(".")
                    prop_key = "properties/"+obj.name+"/"+arr[0]+"/"+arr[1]
                    name = "_".join(prop_key.split("/")[1:])
                    shape = f[prop_key].shape
                    dimensions = self.createDimension(shape)
                    hdf_loc = network.name + ".h5:"+prop_key
                    scalar = self.create_data_item(value=hdf_loc,
                                                   Dimensions=str(dimensions),
                                                   Format='HDF', Precision='8',
                                                   NumberType='Float')
                    el_attr = self.create_attribute(Name=name)
                    el_attr.append(scalar)
                    grid.append(el_attr)

        grid.append(topo)
        grid.append(geo)
        domain.append(grid)
        root.append(domain)

        with open(filename+'.xdmf', 'w') as f:
            f.write(str(self._START))
            f.write(ET.tostring(root).decode("utf-8"))

    def create_domain(self):
        return ET.Element("Domain")

    def create_geometry(self, GeometryType, **attribs):
        element = ET.Element('Geometry')
        element.attrib.update({'GeometryType': GeometryType})
        element.attrib.update(attribs)
        return element

    def create_topology(self, TopologyType, NumberOfElements, **attribs):
        element = ET.Element('Topology')
        element.attrib.update({'TopologyType': TopologyType,
                               'NumberOfElements': NumberOfElements})
        element.attrib.update(attribs)
        return element

    def create_attribute(self, Name, **attribs):
        element = ET.Element('Attribute')
        element.attrib.update({'Name': Name,
                               'AttributeType': 'Scalar',
                               'Center': 'Node'})
        element.attrib.update(attribs)
        return element

    def create_time(self, type='Single', Value=None):
        element = ET.Element('Attribute')
        if type == 'Single' and Value:
            element.attrib['Value'] = Value
        return element

    def create_grid(self, Name, GridType, **attribs):
        element = ET.Element('Grid')
        element.attrib.update({'Name': Name,
                               'GridType': GridType,
                               'Section': None})
        element.attrib.update(attribs)
        if element.attrib['GridType'] is not 'Subset':
            if 'Section' in element.attrib.keys():
                del element.attrib['Section']
        return element

    def create_data_item(self, value, Dimensions, **attribs):
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
