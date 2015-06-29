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
                                           pores=[0])
        except:
            flag = True
        assert flag

    def test_specify_overlapping_geometry(self):
        flag = False
        try:
            OpenPNM.Physics.GenericPhysics(network=self.net,
                                           geometry=self.geo1)
        except:
            flag = True
        assert flag

    def test_get_item_self_name(self):
        a = self.phys1.get('pore.'+self.phys1.name)
        assert a is None
        a = self.phys1['pore.'+self.name]
        assert sp.sum(a) == self.phys1.Np


