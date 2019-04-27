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


if __name__ == '__main__':

    t = GenericPhysicsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
