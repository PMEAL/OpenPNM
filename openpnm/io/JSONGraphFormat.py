import os
import json
import pickle
import numpy as np
from pathlib import Path
from openpnm.utils import logging
from openpnm.io import GenericIO
from openpnm.geometry import Imported
import openpnm.models.geometry as gmods
from openpnm.network import GenericNetwork
logger = logging.getLogger(__name__)


class JSONGraphFormat(GenericIO):
    r"""
    Class for reading and writing OpenPNM networks to JSON Graph Format (JGF).

    Notes
    -----
    The JGF standard must contain data formatted according to
    http://jsongraphformat.info and enforced by JSON schema validation.

    Users must transfer any phase data to the network manually if they wish to
    keep it.
    """

    @classmethod
    def __validate_json__(self, json_file):
        import jsonschema
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
    def save(cls, *args, **kwargs):
        r"""
        This method will be deprecated. Use ``export_data`` instead.
        """
        cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, network, filename=''):
        r"""
        Write the network to disk as a JGF file.

        Parameters
        ----------
        network : OpenPNM Network Object

        filename : string
            Desired file name, defaults to network name if not given
        """

        # Ensure output file is valid
        filename = cls._parse_filename(filename=filename, ext='json')

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

    @classmethod
    def load(cls, *args, **kwargs):
        r"""
        This method will be deprecated.  Use ``import_data`` instead
        """
        return cls.import_data(*args, **kwargs)

    @classmethod
    def import_data(cls, filename, project=None):
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
        filename = cls._parse_filename(filename=filename, ext='json')

        # Load and validate input JSON
        with open(filename, 'r') as file:
            json_file = json.load(file)
            if not cls.__validate_json__(json_file):
                raise Exception('FIle is not in the JSON Graph Format')

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
        network = GenericNetwork(Np=number_of_nodes, Nt=number_of_links)

        # Define primitive throat properties
        network['throat.length'] = link_length
        network['throat.conns'] = np.column_stack([source, target])
        network['throat.diameter'] = 2.0 * np.sqrt(link_squared_radius)

        # Define primitive pore properties
        network['pore.index'] = np.arange(number_of_nodes)
        network['pore.coords'] = np.column_stack([x, y, z])
        network['pore.diameter'] = np.zeros(number_of_nodes)

        geom = Imported(network=network)

        # Define derived throat properties
        geom.add_model(propname='throat.area',
                       model=gmods.throat_cross_sectional_area.cylinder)
        geom.add_model(propname='throat.volume',
                       model=gmods.throat_volume.cylinder)
        geom.add_model(propname='throat.perimeter',
                       model=gmods.throat_perimeter.cylinder)
        geom.add_model(propname='throat.surface_area',
                       model=gmods.throat_surface_area.cylinder)

        # Define derived pore properties
        geom.add_model(propname='pore.area',
                       model=gmods.pore_cross_sectional_area.sphere)
        geom.add_model(propname='pore.volume',
                       model=gmods.pore_volume.sphere)

        return network.project
