import json
import os
import pickle
from pathlib import Path

import jsonschema
import scipy as sp

from openpnm.utils import logging
from openpnm.io import GenericIO
from openpnm.models.geometry import (pore_area, pore_volume, throat_area,
                                     throat_perimeter, throat_surface_area,
                                     throat_volume)
from openpnm.network import GenericNetwork

logger = logging.getLogger(__name__)


class JSONGraphFormat(GenericIO):
    r"""
    Class for reading and writing OpenPNM networks to JSON Graph Format (JGF).

    Notes
    -----
    The JGF standard must contain data formatted according to
    http://jsongraphformat.info and enforced by JSON schema validation.

    Users must transfer any phase data to the network manually if they wish to keep it.
    """

    @classmethod
    def __validate_json__(self, json_file):
        # Validate name of schema file
        relative_path = '../../utils/jgf_schema.pkl'
        schema_file = Path(os.path.realpath(__file__), relative_path)
        schema_file = self._parse_filename(filename=schema_file, ext='pkl')

        # Load schema from pickle file
        with open(schema_file, 'rb') as file:
            jgf_schema = pickle.load(file)

        # Validate JSON agains schema
        try:
            jsonschema.validate(json_file, jgf_schema)
            return True
        except jsonschema.exceptions.ValidationError:
            return False

    @classmethod
    def save(self, network, filename=''):
        r"""
        Write the wetwork to disk as a JGF file.

        Parameters
        ----------
        network : OpenPNM Network Object

        filename : string
            Desired file name, defaults to network name if not given
        """

        # Ensure output file is valid
        filename = self._parse_filename(filename=filename, ext='json')

        # Ensure network contains the required properties
        try:
            required_props = {'pore.diameter', 'pore.coords', 'throat.length',
                              'throat.conns', 'throat.diameter'}
            assert required_props.issubset(network.props())
        except AssertionError:
            raise Exception('Error - network is missing one of: ' +
                            str(required_props))

        # Create 'metadata' JSON object
        graph_metadata_obj = {'number_of_nodes': network.Np,
                              'number_of_links': network.Nt}

        # Create 'nodes' JSON object
        nodes_obj = [
            {
                'id': str(ps),
                'metadata': {
                    'node_squared_radius': int(network['pore.diameter'][ps]/2)**2,
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
                    'link_squared_radius': float(network['throat.diameter'][ts]/2)**2
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

    @classmethod
    def load(self, filename, project=None):
        r"""
        Loads the JGF file onto the given project.

        Parameters
        ----------
        filename : string
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project object is supplied then one will be created and
            returned.

        Returns
        -------
        If no project object is supplied then one will be created and returned.
        """

        # Ensure input file is valid
        filename = self._parse_filename(filename=filename, ext='json')

        # Load and validate input JSON
        with open(filename, 'r') as file:
            json_file = json.load(file)
            if not self.__validate_json__(json_file):
                raise Exception(f'Error - {filename} is not in the JSON Graph Format.')

        # Extract graph metadata from JSON
        number_of_nodes = json_file['graph']['metadata']['number_of_nodes']
        number_of_links = json_file['graph']['metadata']['number_of_links']

        # Extract node properties from JSON
        nodes = sorted(json_file['graph']['nodes'], key=lambda node: int(node['id']))
        x = sp.array([node['metadata']['node_coordinates']['x'] for node in nodes])
        y = sp.array([node['metadata']['node_coordinates']['y'] for node in nodes])
        z = sp.array([node['metadata']['node_coordinates']['z'] for node in nodes])

        # Extract link properties from JSON
        edges = sorted(json_file['graph']['edges'], key=lambda edge: int(edge['id']))
        source = sp.array([int(edge['source']) for edge in edges])
        target = sp.array([int(edge['target']) for edge in edges])
        link_length = sp.array([edge['metadata']['link_length'] for edge in edges])
        link_squared_radius = sp.array(
            [edge['metadata']['link_squared_radius'] for edge in edges])

        # Generate network object
        network = GenericNetwork(Np=number_of_nodes, Nt=number_of_links)

        # Define primitive throat properties
        network['throat.length'] = link_length
        network['throat.conns'] = sp.column_stack([source, target])
        network['throat.diameter'] = 2.0 * sp.sqrt(link_squared_radius)

        # Define derived throat properties
        network['throat.area'] = throat_area.cylinder(network)
        network['throat.volume'] = throat_volume.cylinder(network)
        network['throat.perimeter'] = throat_perimeter.cylinder(network)
        network['throat.surface_area'] = throat_surface_area.cylinder(network)

        # Define primitive pore properties
        network['pore.index'] = sp.arange(number_of_nodes)
        network['pore.coords'] = sp.column_stack([x, y, z])
        network['pore.diameter'] = sp.zeros(number_of_nodes)

        # Define derived pore properties
        network['pore.area'] = pore_area.sphere(network)
        network['pore.volume'] = pore_volume.sphere(network)

        return network.project
