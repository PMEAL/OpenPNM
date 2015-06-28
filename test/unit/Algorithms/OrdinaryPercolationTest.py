import OpenPNM
ctrl = OpenPNM.Base.Controller()
ctrl.loglevel = 60


class OrdinaryPercolationTest:
    def setup_test(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.phase = OpenPNM.Phases.Water(network=self.net)
        self.phys = OpenPNM.Physics.Standard(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
