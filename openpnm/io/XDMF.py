from flatdict import FlatDict
import xml.etree.cElementTree as ET
from openpnm.io import Dict, GenericIO
from openpnm.utils import logging
logger = logging.getLogger(__name__)


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
    def save(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``export_data`` instead.
        """
        cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, network, phases=[], filename=''):
        r"""
        Saves (transient/steady-state) data from the given objects into the
        specified file.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)
            A list of phase objects whose data are to be included

        Notes
        -----
        This method only saves the data, not any of the pore-scale models or
        other attributes.

        """
        import h5py
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        network = network[0]
        # Check if any of the phases has time series
        transient = GenericIO._is_transient(phases=phases)

        if filename == '':
            filename = project.name
        path = cls._parse_filename(filename=filename, ext='xmf')
        # Path is a pathlib object, so slice it up as needed
        fname_xdf = path.name
        d = Dict.to_dict(network, phases=phases, interleave=True,
                         flatten=False, categorize_by=['element', 'data'])
        D = FlatDict(d, delimiter='/')
        # Identify time steps
        t_steps = []
        if transient:
            for key in D.keys():
                if '@' in key:
                    t_steps.append(key.split('@')[1])
        t_steps = list(set(t_steps))
        t_grid = create_grid(Name="TimeSeries", GridType="Collection",
                             CollectionType="Temporal")
        # If steady-state, define '0' time step
        if not transient:
            t_steps.append('0')
        # Setup xdmf file
        root = create_root('Xdmf')
        domain = create_domain()
        # Iterate over time steps present
        for i, t_step in enumerate(t_steps):
            # Define the hdf file
            if not transient:
                fname_hdf = path.stem+".hdf"
            else:
                fname_hdf = path.stem+'@'+t_step+".hdf"
            path_p = path.parent
            f = h5py.File(path_p.joinpath(fname_hdf), "w")
            # Add coordinate and connection information to top of HDF5 file
            f["coordinates"] = network["pore.coords"]
            f["connections"] = network["throat.conns"]
            # geometry coordinates
            row, col = f["coordinates"].shape
            dims = ' '.join((str(row), str(col)))
            hdf_loc = fname_hdf + ":coordinates"
            geo_data = create_data_item(value=hdf_loc, Dimensions=dims,
                                        Format='HDF', DataType="Float")
            geo = create_geometry(GeometryType="XYZ")
            geo.append(geo_data)
            # topolgy connections
            row, col = f["connections"].shape  # col first then row
            dims = ' '.join((str(row), str(col)))
            hdf_loc = fname_hdf + ":connections"
            topo_data = create_data_item(value=hdf_loc, Dimensions=dims,
                                         Format="HDF", NumberType="Int")
            topo = create_topology(TopologyType="Polyline",
                                   NodesPerElement=str(2),
                                   NumberOfElements=str(row))
            topo.append(topo_data)
            # Make HDF5 file with all datasets, and no groups
            for item in D.keys():
                if D[item].dtype == 'O':
                    logger.warning(item + ' has dtype object,'
                                   + ' will not write to file')
                    del D[item]
                elif 'U' in str(D[item][0].dtype):
                    pass
                elif ('@' in item and t_step == item.split('@')[1]):
                    f.create_dataset(name='/'+item.split('@')[0]+'@t',
                                     shape=D[item].shape,
                                     dtype=D[item].dtype,
                                     data=D[item],
                                     compression="gzip")
                elif ('@' not in item and i == 0):
                    f.create_dataset(name='/'+item, shape=D[item].shape,
                                     dtype=D[item].dtype, data=D[item],
                                     compression="gzip")
            # Create a grid
            grid = create_grid(Name=t_step, GridType="Uniform")
            time = create_time(mode='Single', Value=t_step)
            grid.append(time)
            # Add pore and throat properties
            for item in D.keys():
                if item not in ['coordinates', 'connections']:
                    if ("@" in item and t_step == item.split("@")[1]) or (
                        "@" not in item
                    ):
                        attr_type = 'Scalar'
                        shape = D[item].shape
                        dims = (''.join([str(i) + ' ' for i in list(shape)[::-1]]))
                        if '@' in item:
                            item = item.split('@')[0]+'@t'
                            hdf_loc = fname_hdf + ":" + item
                        elif ('@' not in item and i == 0):
                            hdf_loc = fname_hdf + ":" + item
                        elif ('@' not in item and i > 0):
                            hdf_loc = path.stem + '@' + t_steps[0] + ".hdf" + ":" + item
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
                    else:
                        pass
            grid.append(topo)
            grid.append(geo)
            t_grid.append(grid)
            # CLose the HDF5 file
            f.close()
        domain.append(t_grid)
        root.append(domain)
        with open(path_p.joinpath(fname_xdf), 'w') as file:
            file.write(cls._header)
            file.write(ET.tostring(root).decode("utf-8"))


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


def create_time(mode='Single', Value=None):
    element = ET.Element('Time')
    if mode == 'Single' and Value:
        element.attrib['Value'] = Value
    return element


def create_grid(Name, GridType, **attribs):
    element = ET.Element('Grid')
    element.attrib.update({'Name': Name,
                           'GridType': GridType,
                           'Section': None})
    element.attrib.update(attribs)
    if element.attrib['GridType'] != 'Subset':
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
