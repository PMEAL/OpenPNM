import openpnm as op
from numpy.testing import assert_allclose
import pytest


class PhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_instantiate_without_network_fails(self):
        with pytest.raises(TypeError):
            op.phase.Phase()

    def test_getitem_interpolation(self):
        pn = op.network.Demo()
        air = op.phase.Air(network=pn)
        assert_allclose(air['pore.density'][0], 1.17982341)
        assert_allclose(air['throat.density'][0], 1.17982341)
        air['throat.test'] = 0.5
        assert_allclose(air['pore.test'][0], 0.5)
        with pytest.raises(KeyError):
            air['throat.density1']
        with pytest.raises(KeyError):
            air['pore.density1']

    def test_getitem_from_network(self):
        pn = op.network.Demo()
        air = op.phase.Air(network=pn)
        assert air['pore.left'].sum() == 3

    def test_getitem_domains(self):
        pn = op.network.Demo()
        air = op.phase.Air(network=pn)
        assert len(air['pore.density@left']) == 3
        assert_allclose(air['pore.density@left'][0], 1.17982341)
        air['pore.test'] = False
        air['pore.test'][[0, 1]] = True
        assert len(air['pore.density@test']) == 2
        with pytest.raises(KeyError):
            air['pore.density@blah']

    def test_getitem_param(self):
        pn = op.network.Demo()
        air = op.phase.Air(network=pn)
        assert len(air['param.formula']) == 2

    def test_getitem_interpolation_disabled(self):
        pn = op.network.Demo()
        air = op.phase.Air(network=pn)
        air.settings['auto_interpolate'] = False
        with pytest.raises(KeyError):
            air['throat.density']


if __name__ == '__main__':

    t = PhaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
