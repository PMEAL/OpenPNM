import os
import json
import pickle
import logging
import numpy as np
from pathlib import Path
from openpnm.io import _parse_filename
from openpnm.network import Network
logger = logging.getLogger(__name__)


def _validate_json(json_file):
    import jsonschema
    # Validate name of schema file
    relative_path = '../../utils/jgf_schema.pkl'
    schema_file = Path(os.path.realpath(__file__), relative_path)
    schema_file = _parse_filename(filename=schema_file, ext='pkl')

    # Load schema from pickle file
    with open(schema_file, 'rb') as file:
        jgf_schema = pickle.load(file)

    # Validate JSON agains schema
    try:
        jsonschema.validate(json_file, jgf_schema)
        return True
    except jsonschema.exceptions.ValidationError:
        return False


def network_to_jsongraph(network, filename=''):
    r"""
    Write the network to disk as a JGF file.

    Parameters
    ----------
    network : Network

    filename : str
        Desired file name, defaults to network name if not given
    """

    # Ensure output file is valid
    filename = _parse_filename(filename=filename, ext='json')

    # Ensure network contains the required properties
    try:
        required_props = {'pore.diameter', 'pore.coords', 'throat.length',
                          'throat.conns', 'throat.diameter'}
        assert required_props.issubset(network.props())
    except AssertionError:
        raise Exception('Error - network is missing one of: '
                        + str(required_props))

    # Create 'metadata' JSON object
    graph_metadata_obj = {'number_of_nodes': network.Np,
                          'number_of_links': network.Nt}

    # Create 'nodes' JSON object
    nodes_obj = [
        {
            'id': str(ps),
            'metadata': {
                'node_squared_radius': int(network['pore.diameter'][ps] / 2)**2,
                'node_coordinates': {
                    'x': int(network['pore.coords'][ps, 0]),
                    'y': int(network['pore.coords'][ps, 1]),
                    'z': int(network['pore.coords'][ps, 2])
                }
            }
        } for ps in network.Ps]

    # Create 'edges' JSON object
    edges_obj = [
        {
            'id': str(ts),
            'source': str(network['throat.conns'][ts, 0]),
            'target': str(network['throat.conns'][ts, 1]),
            'metadata': {
                'link_length': float(network['throat.length'][ts]),
                'link_squared_radius': float(network['throat.diameter'][ts] / 2)**2
            }
        } for ts in network.Ts]

    # Build 'graph' JSON object from 'metadata', 'nodes' and 'edges'
    graph_obj = {'metadata': graph_metadata_obj,
                 'nodes': nodes_obj,
                 'edges': edges_obj}

    # Build full JSON object
    json_obj = {'graph': graph_obj}

    # Load and validate input JSON
    with open(filename, 'w') as file:
        json.dump(json_obj, file, indent=2)


def network_from_jsongraph(filename):
    r"""
    Loads the JGF file onto the given project.

    Parameters
    ----------
    filename : str
        The name of the file containing the data to import.  The formatting
        of this file is outlined below.

    Returns
    -------
    network : dict
        An OpenPNM Network dictionary

    """

    # Ensure input file is valid
    filename = _parse_filename(filename=filename, ext='json')

    # Load and validate input JSON
    with open(filename, 'r') as file:
        json_file = json.load(file)
        if not _validate_json(json_file):
            raise Exception('File is not in the JSON Graph Format')

    # Extract graph metadata from JSON
    number_of_nodes = json_file['graph']['metadata']['number_of_nodes']
    number_of_links = json_file['graph']['metadata']['number_of_links']

    # Extract node properties from JSON
    nodes = sorted(json_file['graph']['nodes'], key=lambda node: int(node['id']))
    x = np.array([node['metadata']['node_coordinates']['x'] for node in nodes])
    y = np.array([node['metadata']['node_coordinates']['y'] for node in nodes])
    z = np.array([node['metadata']['node_coordinates']['z'] for node in nodes])

    # Extract link properties from JSON
    edges = sorted(json_file['graph']['edges'], key=lambda edge: int(edge['id']))
    source = np.array([int(edge['source']) for edge in edges])
    target = np.array([int(edge['target']) for edge in edges])
    link_length = np.array([edge['metadata']['link_length'] for edge in edges])
    link_squared_radius = np.array(
        [edge['metadata']['link_squared_radius'] for edge in edges])

    # Generate network object
    network = Network()
    network['pore.all'] = np.ones([number_of_nodes, ], dtype=bool)
    network['throat.all'] = np.ones([number_of_links, ], dtype=bool)

    # Define primitive throat properties
    network['throat.length'] = link_length
    network['throat.conns'] = np.column_stack([source, target])
    network['throat.diameter'] = 2.0 * np.sqrt(link_squared_radius)

    # Define primitive pore properties
    network['pore.index'] = np.arange(number_of_nodes)
    network['pore.coords'] = np.column_stack([x, y, z])
    network['pore.diameter'] = np.zeros(number_of_nodes)

    return network
