import OpenPNM


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
