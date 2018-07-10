import copy
import os
from pathlib import Path

import py
import pytest
import scipy as sp

import openpnm as op


class JSONGraphFormatTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net['pore.diameter'] = 2.0 * sp.ones(self.net.Np)
        self.net['throat.diameter'] = 2.0 * sp.ones(self.net.Nt)
        self.net.add_model(propname='throat.length',
                           model=op.models.geometry.throat_length.ctc)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_validation_success(self):
        json_obj = {'graph': {'nodes': [{'id': "0"}]}}  # 'id' is a string
        jgf = op.io.JSONGraphFormat()
        assert jgf.__validate_json__(json_obj)

    def test_validation_failure(self):
        json_obj = {'graph': {'nodes': [{'id': 0}]}}    # 'id' is not a string
        jgf = op.io.JSONGraphFormat()
        assert not jgf.__validate_json__(json_obj)

    def test_save_failure(self):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/JSONGraphFormat')
        filename = Path(path.resolve(), 'save_failure.json')

        # Create a deep copy of network with one required property missing
        net = copy.deepcopy(self.net)
        net.pop('pore.diameter')

        # Ensure an exception was thrown
        with pytest.raises(Exception) as e_info:
            op.io.JSONGraphFormat.save(net, filename=filename)
        expected_error = 'Error - network is missing one of:'
        assert expected_error in str(e_info.value)

    def test_load_failure(self):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/JSONGraphFormat')
        filename = Path(path.resolve(), 'invalid.json')

        # Ensure an exception was thrown
        with pytest.raises(Exception) as e_info:
            op.io.JSONGraphFormat.load(filename)
        expected_error = f'Error - {filename} is not in the JSON Graph Format.'
        assert expected_error in str(e_info.value)

    def test_load_success(self):
        # Load JSON file and ensure project integrity
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/JSONGraphFormat')
        filename = Path(path.resolve(), '2nodes_1link.json')
        project = op.io.JSONGraphFormat.load(filename)
        assert len(project) == 1

        # Ensure overal network properties
        net = project.network
        assert net.Np == 2
        assert net.Nt == 1

        # Ensure existence of pore properties
        pore_props = {'pore.index', 'pore.coords', 'pore.diameter',
                      'pore.area', 'pore.volume'}
        assert pore_props.issubset(net.props())

        # Ensure correctness of pore properties
        assert sp.array_equal(net['pore.area'], sp.array([0, 0]))
        assert sp.array_equal(net['pore.index'], sp.array([0, 1]))
        assert sp.array_equal(net['pore.volume'], sp.array([0, 0]))
        assert sp.array_equal(net['pore.diameter'], sp.array([0, 0]))
        assert sp.array_equal(net['pore.coords'][0], sp.array([0, 0, 0]))
        assert sp.array_equal(net['pore.coords'][1], sp.array([1, 1, 1]))

        # Ensure existence of throat properties
        throat_props = {'throat.length', 'throat.conns', 'throat.diameter',
                        'throat.area', 'throat.volume', 'throat.perimeter',
                        'throat.surface_area'}
        assert throat_props.issubset(net.props())

        # Ensure correctness of throat properties
        length = 1.73205080757
        squared_radius = 5.169298742047715
        assert net['throat.length'] == length
        assert net['throat.area'] == sp.pi * squared_radius
        assert sp.array_equal(net['throat.conns'], sp.array([[0, 1]]))
        assert net['throat.diameter'] == 2.0 * sp.sqrt(squared_radius)
        assert net['throat.volume'] == sp.pi * squared_radius * length
        assert net['throat.perimeter'] == 2.0 * sp.pi * sp.sqrt(squared_radius)
        assert net['throat.surface_area'] == 2.0 * \
            sp.sqrt(squared_radius) * sp.pi * length


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = JSONGraphFormatTest()
    self = t  # For interacting with the tests at the command line
    tmpdir = py.path.local()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=tmpdir)
