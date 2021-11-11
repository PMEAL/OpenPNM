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
        assert 'throat.conduit_lengths' not in geo.keys()
        assert 'throat.length' in geo.models.keys()
        assert 'throat.diffusive_size_factors' in geo.models.keys()
        assert 'throat.hydraulic_size_factors' in geo.models.keys()

    def test_with_added_pores(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        net['pore.diameter'] = 2.0
        net['throat.diameter'] = 1.0
        _ = op.geometry.Imported(network=net)
        op.topotools.extend(network=net,
                            pore_coords=[[1, 1, 1]],
                            throat_conns=[[0, 27]])
        h = net.project.check_geometry_health()
        assert h.health is False
        _ = op.geometry.GenericGeometry(network=net, pores=27, throats=54)
        h = net.project.check_geometry_health()
        assert h.health is True
        assert np.any(np.isnan(net['pore.diameter']))
        assert np.any(np.isnan(net['throat.diameter']))

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
