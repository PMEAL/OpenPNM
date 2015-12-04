import OpenPNM
import scipy as sp


class GenericPhysicsTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        P1 = self.net.pores('top', mode='not')
        T1 = self.net.find_neighbor_throats(pores=P1)
        self.geo1 = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                     pores=P1,
                                                     throats=T1)
        P2 = self.net.pores('top')
        T2 = self.net.find_neighbor_throats(pores=P2, mode='intersection')
        self.geo2 = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                     pores=P2,
                                                     throats=T2)
        self.phase1 = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase2 = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phys1 = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                    phase=self.phase1,
                                                    geometry=self.geo1)

    def test_specify_pores_and_geometry(self):
        flag = False
        try:
            OpenPNM.Physics.GenericPhysics(network=self.net,
                                           phase=self.phase1,
                                           geometry=self.geo1,
                                           pores=[0])
        except:
            flag = True
        assert flag

    def test_specify_overlapping_pores(self):
        flag = False
        try:
            OpenPNM.Physics.GenericPhysics(network=self.net,
                                           phase=self.phase1,
                                           pores=[0])
        except:
            flag = True
        assert flag

    def test_specify_overlapping_geometry(self):
        flag = False
        try:
            OpenPNM.Physics.GenericPhysics(network=self.net,
                                           phase=self.phase1,
                                           geometry=self.geo1)
        except:
            flag = True
        assert flag

    def test_get_item_self_name(self):
        a = self.phys1.get('pore.'+self.phys1.name)
        assert a is None
        a = self.phys1['pore.'+self.phys1.name]
        assert sp.sum(a) == self.phys1.Np

    def test_instantiation_with_no_geom_or_phase(self):
        phys = OpenPNM.Physics.GenericPhysics(network=self.net)
        ctrl = phys.controller
        assert phys in ctrl.values()
        phase = phys._phases[0]
        assert phase in ctrl.values()
        assert phys.name in phase.physics()
        assert phase.name in phys.phases()

    def test_reassign_to_different_phase(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        geo1 = OpenPNM.Geometry.GenericGeometry(network=net,
                                                pores=net.Ps,
                                                throats=net.Ts)
        phase1 = OpenPNM.Phases.GenericPhase(network=net)
        phase2 = OpenPNM.Phases.GenericPhase(network=net)
        phys = OpenPNM.Physics.GenericPhysics(network=net,
                                              phase=phase1,
                                              geometry=geo1)

        assert phys.parent_phase is phase1
        phys.parent_phase = phase2
        assert 'pore.'+phys.name not in phase1.labels()
        assert 'pore.'+phys.name in phase2.labels()
        assert phys.phases()[0] == phase2.name
        assert phase2.physics()[0] == phys.name
        a = phase2.check_physics_health()
        assert sum([len(value) for value in a.values()]) == 0
        a = phase1.check_physics_health()
        assert len(a['undefined_pores']) > 0
