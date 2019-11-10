import openpnm as op
import pytest
ws = op.Workspace()
ws.settings['loglevel'] = 10


class GenericPhysicsTest:

    def setup_class(self):
        pass

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_instantiate_normal(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                         geometry=geo)
        assert [net, geo, phase, phys] == phase.project

    def test_instantiate_with_only_network(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        phys = op.physics.GenericPhysics(network=net)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is not None
        assert phys.project.find_geometry(phys) is not None

    def test_instantiate_with_only_network_and_phase(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is phase
        assert phys.project.find_geometry(phys) is not None

    def test_instantiate_with_network_and_geometry(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net)
        phase = op.phases.GenericPhase(network=net)  # Don't use this phase
        phys = op.physics.GenericPhysics(network=net, geometry=geo)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is not phase
        phase = phys.project.find_full_domain(phys)
        assert phys.project.find_phase(phys) is phase
        assert phys.project.find_geometry(phys) is geo

    def test_instantiate_with_physics_present(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net)
        phys = op.physics.GenericPhysics(network=net, geometry=geo)
        with pytest.raises(Exception):
            phys = op.physics.GenericPhysics(network=net)

    def test_instantiate_with_geometry_present(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net)
        phys = op.physics.GenericPhysics(network=net)
        with pytest.raises(Exception):
            phys = op.physics.GenericPhysics(network=net)

    def test_set_phase(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net,
                                          pores=net.Ps, throats=net.Ts)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                         geometry=geo)
        phase2 = op.phases.GenericPhase(network=net)
        with pytest.raises(Exception):
            phys.set_phase(phase=phase2, mode='add')
        with pytest.raises(Exception):
            phys.set_phase(phase=phase2, mode='remove')
        phys.set_phase(phase=phase, mode='remove')
        phys.set_phase(phase=phase2, mode='add')
        phys.set_phase(phase=phase, mode='swap')

    def test_set_geometry(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        geo1 = op.geometry.GenericGeometry(network=net, pores=[0, 1, 2, 3],
                                           throats=net.Ts)
        geo2 = op.geometry.GenericGeometry(network=net, pores=[4, 5, 6, 7])
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                         geometry=geo1)
        phys.set_geometry(geometry=geo2)

    def test_set_phase_and_geometry_from_different_project(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        net2 = op.network.Cubic(shape=[5, 5, 5])
        geo2 = op.geometry.GenericGeometry(network=net2)
        phase2 = op.phases.GenericPhase(network=net2)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net)
        with pytest.raises(Exception):
            phys.set_phase(phase=phase2)
        phys.set_phase(phase=phase, mode='swap')
        with pytest.raises(Exception):
            phys.set_geometry(geometry=geo2)


if __name__ == '__main__':

    t = GenericPhysicsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
