import openpnm as op
import numpy as np
import pytest
ws = op.Workspace()
ws.settings['loglevel'] = 50


class ImportedTest:

    def setup_class(self):
        pass

    def test_init(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        net['pore.diameter'] = 2.0
        net['throat.diameter'] = 1.0
        geo = op.geometry.Imported(network=net)
        assert 'pore.diameter' in geo.keys()
        assert 'pore.diameter' not in net.keys()
        assert 'throat.diameter' in geo.keys()
        assert 'throat.diameter' not in net.keys()

    def test_init_no_props_on_network(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.Imported(network=net)
        assert 'throat.length' not in geo.keys()
        assert 'throat.endpoints' not in geo.keys()
        assert 'throat.conduit_lengths' not in geo.keys()
        assert 'throat.length' in geo.models.keys()
        assert 'throat.endpoints' in geo.models.keys()
        assert 'throat.conduit_lengths' in geo.models.keys()

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()


if __name__ == '__main__':

    t = ImportedTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
