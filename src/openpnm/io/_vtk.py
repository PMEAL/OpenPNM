import logging
import numpy as np
import pandas as pd
from xml.etree import ElementTree as ET
from openpnm.io import project_to_dict, _parse_filename
from openpnm.utils import Workspace
from openpnm.utils._misc import is_transient
logger = logging.getLogger(__name__)
ws = Workspace()


_TEMPLATE = """
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
""".strip()


def project_to_vtk(project, filename="",
                   fill_nans=None, fill_infs=None):
    r"""
    Save network and phase data to a single vtp file for visualizing in
    Paraview.

    Parameters
    ----------
    network : Network
        The Network containing the data to be written
    phases : list, optional
        A list containing Phase(s) containing data to be
        written
    filename : str, optional
        Filename to write data.  If no name is given the file is named
        after the network
    fill_nans : scalar
        The value to use to replace NaNs with.  The VTK file format does
        not work with NaNs, so they must be dealt with.  The default is
        `None` which means property arrays with NaNs are not written to the
        file.  Other useful options might be 0 or -1, but the user must
        be aware that these are not real values, only place holders.
    fill_infs : scalar
        The value to use to replace infs with.  The default is ``None``
        which means that property arrays containing ``None`` will *not*
        be written to the file, and a warning will be issued.  A useful
        value is

    """
    network = project.network
    algs = project.algorithms
    # Check if any of the phases has time series
    transient = is_transient(algs)
    if transient:
        logger.warning(
            "vtp format does not support transient data, " + "use xdmf instead"
        )
    if filename == "":
        filename = project.name
    filename = _parse_filename(filename=filename, ext="vtp")

    am = project_to_dict(project=project,
                         categorize_by=["object", "data"])
    am = pd.json_normalize(am, sep='.').to_dict(orient='records')[0]
    for k in list(am.keys()):
        am[k.replace('.', ' | ')] = am.pop(k)
    key_list = list(sorted(am.keys()))

    points = network["pore.coords"]
    pairs = network["throat.conns"]
    num_points = np.shape(points)[0]
    num_throats = np.shape(pairs)[0]

    root = ET.fromstring(_TEMPLATE)
    piece_node = root.find("PolyData").find("Piece")
    piece_node.set("NumberOfPoints", str(num_points))
    piece_node.set("NumberOfLines", str(num_throats))
    points_node = piece_node.find("Points")
    coords = _array_to_element("coords", points.T.ravel("F"), n=3)
    points_node.append(coords)
    lines_node = piece_node.find("Lines")
    connectivity = _array_to_element("connectivity", pairs)
    lines_node.append(connectivity)
    offsets = _array_to_element("offsets", 2 * np.arange(len(pairs)) + 2)
    lines_node.append(offsets)

    point_data_node = piece_node.find("PointData")
    cell_data_node = piece_node.find("CellData")
    for key in key_list:
        array = am[key]
        if array.dtype == "O":
            logger.warning(key + " has dtype object," + " will not write to file")
        else:
            if array.dtype == bool:
                array = array.astype(int)
            if np.any(np.isnan(array)):
                if fill_nans is None:
                    logger.warning(key + " has nans," + " will not write to file")
                    continue
                else:
                    array[np.isnan(array)] = fill_nans
            if np.any(np.isinf(array)):
                if fill_infs is None:
                    logger.warning(key + " has infs," + " will not write to file")
                    continue
                else:
                    array[np.isinf(array)] = fill_infs
            element = _array_to_element(key, array)
            if array.size == num_points:
                point_data_node.append(element)
            elif array.size == num_throats:
                cell_data_node.append(element)

    tree = ET.ElementTree(root)
    tree.write(filename)

    with open(filename, "r+") as f:
        string = f.read()
        string = string.replace("</DataArray>", "</DataArray>\n\t\t\t")
        f.seek(0)
        # consider adding header: '<?xml version="1.0"?>\n'+
        f.write(string)


def _array_to_element(name, array, n=1):
    dtype_map = {
        "int8": "Int8",
        "int16": "Int16",
        "int32": "Int32",
        "int64": "Int64",
        "uint8": "UInt8",
        "uint16": "UInt16",
        "uint32": "UInt32",
        "uint64": "UInt64",
        "float32": "Float32",
        "float64": "Float64",
        "str": "String",
    }
    element = None
    if str(array.dtype) in dtype_map.keys():
        element = ET.Element("DataArray")
        element.set("Name", name)
        element.set("NumberOfComponents", str(n))
        element.set("type", dtype_map[str(array.dtype)])
        element.text = "\t".join(map(str, array.ravel()))
    return element


def _element_to_array(element, n=1):
    string = element.text
    dtype = element.get("type")
    array = np.fromstring(string, sep="\t")
    array = array.astype(dtype.lower())
    if n != 1:
        array = array.reshape(array.size // n, n)
    return array
