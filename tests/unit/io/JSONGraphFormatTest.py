import os
from pathlib import Path

import openpnm as op
import py


class JSONGraphFormatTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_load(self):
        # Load JSON file and ensure project integrity
        path = Path(os.path.realpath(__file__), '../../fixtures/JSONGraphFormat')
        project = op.io.JSONGraphFormat.load(path + '/2nodes_1link.json')
        assert len(project) == 1

        # Ensure overal network properties
        net = project.network
        assert net.Np == 2
        assert net.Nt == 1

        # Ensure existence of throat properties
        throat_props = {'throat.length', 'throat.conns', 'throat.diameter',
                        'throat.area', 'throat.volume', 'throat.perimeter',
                        'throat.surface_area'}
        assert throat_props.issubset(net.props())

        # Ensure existence of pore properties
        pore_props = {'pore.index', 'pore.coords', 'pore.diameter',
                      'pore.area', 'pore.volume'}
        assert pore_props.issubset(net.props())


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
