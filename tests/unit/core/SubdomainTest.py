import pytest
import numpy as np
import openpnm as op


class SubdomainTest:

    def setup_class(self):
        pass

    def test_add_locations(self):
        pn = op.network.Cubic([6, 1, 1])
        g1 = op.geometry.GenericGeometry(network=pn)
        g2 = op.geometry.GenericGeometry(network=pn)
        g1.set_locations(pores=[0, 1, 2], mode='add')
        g2.set_locations(pores=[3, 4, 5], mode='add')
        assert np.all(g1.Ps == [0, 1, 2])
        assert np.all(g2.Ps == [0, 1, 2])

    def test_fail_to_switch_locations(self):
        pn = op.network.Cubic([6, 1, 1])
        g1 = op.geometry.GenericGeometry(network=pn)
        g2 = op.geometry.GenericGeometry(network=pn)
        g1.set_locations(pores=[0, 1, 2], mode='add')
        with pytest.raises(Exception):
            g2.set_locations(pores=[0, 1, 2], mode='add')
        assert np.all(g1.Ps == [0, 1, 2])
        assert np.all(g2.Ps == [ ])

    def test_switch_locations(self):
        pn = op.network.Cubic([6, 1, 1])
        g1 = op.geometry.GenericGeometry(network=pn)
        g2 = op.geometry.GenericGeometry(network=pn)
        g1.set_locations(pores=[0, 1, 2], mode='add')
        assert np.all(g1.Ps == [0, 1, 2])
        assert np.all(g2.Ps == [ ])
        g2.set_locations(pores=[0, 1, 2], mode='switch')
        assert np.all(g1.Ps == [ ])
        assert np.all(g2.Ps == [0, 1, 2])

    def test_drop_locations(self):
        pn = op.network.Cubic([6, 1, 1])
        g1 = op.geometry.GenericGeometry(network=pn)
        g2 = op.geometry.GenericGeometry(network=pn)
        g1.set_locations(pores=[0, 1, 2], mode='add')
        assert np.all(g1.Ps == [0, 1, 2])
        assert np.all(g2.Ps == [ ])
        g1.set_locations(pores=[0, 1, 2], mode='drop')
        assert np.all(g1.Ps == [ ])
        assert np.all(g2.Ps == [ ])


if __name__ == '__main__':

    t = SubdomainTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
