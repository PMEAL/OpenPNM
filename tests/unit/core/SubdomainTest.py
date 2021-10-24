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

    def test_drop_pores_from_geo_add_phys(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                          geometry=geo)
        geo.set_locations(pores=[0], mode='drop')
        phys.set_locations(pores=[0], mode='drop')
        geo2 = op.geometry.GenericGeometry(network=net, pores=[0])
        phys2 = op.physics.GenericPhysics(network=net, phase=phase,
                                          geometry=geo2)
        assert geo.Np == 26
        assert phys.Np == 26
        assert geo2.Np == 1
        assert phys2.Np == 1

    def test_phase_and_geom_attr_with_two_domains(self):
        pn = op.network.Cubic([6, 1, 1])
        g1 = op.geometry.GenericGeometry(network=pn, pores=[0, 1, 2])
        g2 = op.geometry.GenericGeometry(network=pn, pores=[3, 4, 5])
        air = op.phases.Air(network=pn)
        phys1 = op.physics.GenericPhysics(network=pn)
        phys2 = op.physics.GenericPhysics(network=pn)
        with pytest.raises(Exception):
            phys1.geometry = g1  # Can't assign a geo if no phase assigned
        phys1.phase = air
        phys1.geometry = g1
        assert phys1 in g1.physics
        phys2.phase = air
        with pytest.raises(Exception):
            phys2.geometry = g1  # g1 is already assigned to phys1
        phys2.geometry = g2
        assert phys2 in g2.physics
        assert phys1 in air.physics
        assert phys2 in air.physics


if __name__ == '__main__':

    t = SubdomainTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
