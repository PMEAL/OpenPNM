import openpnm as op
import pytest


class GenericPhysicsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
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

    def test_instantiate_with_phase_only(self):
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net,
                                         phase=phase)
        assert phys.project is not None
        assert phys.project.find_phase(phys) is phase
        with pytest.raises(Exception):
            phys.project.find_geometry(phys)

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
        phase = op.phases.GenericPhase(network=self.net)
        phys = op.physics.GenericPhysics(network=self.net)
        assert phys.project is not None
        with pytest.raises(Exception):
            phys.project.find_phase(phys) is phase
        with pytest.raises(Exception):
            phys.project.find_geometry(phys)

    def test_set_phase_afer_instantiation(self):
        phys = op.physics.GenericPhysics(network=self.net)
        phase = op.phases.GenericPhase(network=self.net)
        assert 'pore.' + phys.name not in phase.keys()
        phys.set_phase(phase=phase, mode='add')
        assert 'pore.' + phys.name in phase.keys()
        phys.set_phase(phase=phase, mode='remove')
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
