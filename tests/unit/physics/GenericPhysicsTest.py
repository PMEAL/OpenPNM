import openpnm as op
import pytest


class GenericPhysicsTest:

    def setup_class(self):
        pass

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_instantiate_normally(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase,
                                         geometry=geo)

    def test_instantiate_with_phase_only(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, phase=phase)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is phase
        assert phys.project.find_geometry(phys) is not None

    def test_instantiate_with_no_phase(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.GenericGeometry(network=net)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net, geometry=geo)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is not phase
        phase = phys.project.find_full_domain(phys)
        assert phys.project.find_phase(phys) is phase
        assert phys.project.find_geometry(phys) is geo

    def test_instantiate_with_only_network(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        phys = op.physics.GenericPhysics(network=net)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is not None
        assert phys.project.find_geometry(phys) is not None

    def test_set_phase_afer_instantiation(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        phys = op.physics.GenericPhysics(network=net)
        phase = op.phases.GenericPhase(network=net)
        assert 'pore.' + phys.name not in phase.keys()
        phys.set_phase(phase=phase, mode='add')
        assert 'pore.' + phys.name in phase.keys()
        phys.set_phase(phase=phase, mode='remove')
        assert 'pore.' + phys.name not in phase.keys()

    def test_set_geom_after_instantiation(self):
        net = op.network.Cubic(shape=[3, 3, 3])
#        phase = op.phases.GenericPhase(network=net)
#        phys = op.physics.GenericPhysics(network=net, phase=phase)
#        assert phase['pore.'+phys.name].sum() == net.Np
#        phys.set_geometry(geometry=self.geo, pores=[])
#        assert phase['pore.'+phys.name].sum() == 0

    def test_set_phase_and_geometry_from_different_project(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        net2 = op.network.Cubic(shape=[5, 5, 5])
        geo2 = op.geometry.GenericGeometry(network=net2)
        phase2 = op.phases.GenericPhase(network=net2)
        phase = op.phases.GenericPhase(network=net)
        phys = op.physics.GenericPhysics(network=net)
        with pytest.raises(Exception):
            phys.set_phase(phase=phase2)
        phys.set_phase(phase=phase)
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
