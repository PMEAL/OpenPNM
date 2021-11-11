import openpnm as op
import pytest


class GenericPhysicsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_instantiate_normally(self):
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net,
                                         phase=phase,
                                         geometry=self.geo)
        assert phys.Np == self.geo.Np
        assert phys.Nt == self.geo.Nt

    def test_instantiate_with_phase_only(self):
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net,
                                         phase=phase)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is phase
        with pytest.raises(Exception):
            _ = phys.project.find_geometry(phys)

    def test_instantiate_with_geometry_only(self):
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net,
                                         geometry=self.geo)
        assert phys.project is not None
        with pytest.raises(Exception):
            phys.project.find_phase(phys) is phase
        with pytest.raises(Exception):
            phys.project.find_geometry(phys)

    def test_instantiate_with_only_network(self):
        _ = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net)
        assert phys.project is not None
        with pytest.raises(Exception):
            _ = phys.project.find_phase(phys)
        with pytest.raises(Exception):
            _ = phys.project.find_geometry(phys)

    def test_instantiate_with_pores_and_throats(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        _ = op.geometry.GenericGeometry(network=net,
                                        pores=net.Ps,
                                        throats=net.Ts)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net,
                                         phase=phase,
                                         pores=net.Ps,
                                         throats=net.Ts)
        assert phys.Np == 27
        assert phys.Nt == 54

    def test_set_phase_after_instantiation(self):
        phys = op.physics.GenericPhysics(network=self.net)
        phase = op.phases.GenericPhase(network=self.net)
        assert 'pore.' + phys.name not in phase.keys()
        phys.set_phase(phase=phase, mode='add')
        assert 'pore.' + phys.name in phase.keys()
        phys.set_phase(phase=phase, mode='drop')
        assert 'pore.' + phys.name not in phase.keys()

    def test_set_geom_after_instantiation(self):
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net, phase=phase)
        assert phase['pore.'+phys.name].sum() == 0
        phys.set_geometry(geometry=self.geo)
        assert phase['pore.'+phys.name].sum() == phase.Np

    def test_set_phase_and_geometry_from_different_project(self):
        net2 = op.network.Cubic(shape=[5, 5, 5])
        geo2 = op.geometry.GenericGeometry(network=net2)
        phase2 = op.phases.GenericPhase(network=net2)
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net)
        with pytest.raises(Exception):
            _ = phys.set_phase(phase=phase2)
        phys.set_phase(phase=phase)
        with pytest.raises(Exception):
            _ = phys.set_geometry(geometry=geo2)

    def test_swap_phase(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase1 = op.phases.GenericPhase(network=net)
        phase2 = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase1,
                                         geometry=geo)
        phys.set_phase(phase=phase2, mode='move')
        assert phys.Np == 27
        assert phys.Nt == 54

    def test_drop_add_phase(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase1 = op.phases.GenericPhase(network=net)
        phase2 = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase1,
                                         geometry=geo)
        phys.set_phase(mode='drop')
        assert phys.Np == 0
        assert phys.Nt == 0
        phys.set_phase(phase=phase2, mode='add')
        assert phys.Np == 0
        assert phys.Nt == 0

    def test_drop_add_geometry(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                         geometry=geo)
        assert phys.Np == 27
        assert phys.Nt == 54
        phys.set_geometry(mode='drop')
        assert phys.Np == 0
        assert phys.Nt == 0
        phys.set_geometry(geometry=geo, mode='add')
        assert phys.Np == 27
        assert phys.Nt == 54
        phys.set_phase(mode='drop')
        assert phys.Np == 0
        assert phys.Nt == 0
        with pytest.raises(Exception):
            phys.set_geometry(geometry=geo, mode='add')

    def test_using_geometry_attr(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                         geometry=geo)
        assert phys.Np == 27
        del phys.geometry
        assert phys.Np == 0
        phys.geometry = geo
        assert phys.Np == 27
        assert phys.geometry is geo

    def test_using_phase_attr(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase1 = op.phases.GenericPhase(network=net)
        phase2 = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase1,
                                         geometry=geo)
        assert phys.Np == 27
        assert phys in phase1.physics
        assert phys not in phase2.physics
        phys.phase = phase2
        assert phys in phase2.physics
        assert phys not in phase1.physics
        assert phys.Np == 27

    def test_set_phase_modes(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        phase1 = op.phases.GenericPhase(network=net)
        phase2 = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase1,
                                         geometry=geo)
        assert phys.Np == 27
        assert phys.Nt == 54
        assert phys in phase1.physics
        assert phys not in phase2.physics
        with pytest.raises(Exception):
            phys.set_phase(phase=phase2, mode='add')
        phys.set_phase(phase=phase2, mode='move')
        assert phys in phase2.physics
        assert phys not in phase1.physics
        assert phys.Np == 27
        assert phys.Nt == 54
        phys.set_phase(mode='drop')
        phys.set_phase(phase=phase1, mode='add')
        assert phys.Np == 0
        assert phys.Nt == 0
        assert phys in phase1.physics
        assert phys not in phase2.physics

    def test_set_geometry_mode_move(self):
        pn = op.network.Cubic([6, 1, 1])
        g1 = op.geometry.GenericGeometry(network=pn, pores=[0, 1, 2])
        g2 = op.geometry.GenericGeometry(network=pn, pores=[3, 4, 5])
        air = op.phases.Air(network=pn)
        phys1 = op.physics.GenericPhysics(network=pn, phase=air, geometry=g1)
        phys2 = op.physics.GenericPhysics(network=pn, phase=air)
        phys2.set_geometry(geometry=g1, mode='move')
        assert phys1.Np == 0
        assert phys2.Np == 3
        phys1.set_geometry(geometry=g2, mode='add')
        assert phys1.Np == 3
        assert phys2.Np == 3


if __name__ == '__main__':

    t = GenericPhysicsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
