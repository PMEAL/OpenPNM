import OpenPNM
ctrl = OpenPNM.Base.Controller()
ctrl.loglevel = 60


class DrainageTest:
    def setup_test(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.phase = OpenPNM.Phases.Water(network=self.net)
        self.phys = OpenPNM.Physics.Standard(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.alg = OpenPNM.Algorithms.Drainage(network=self.net)

    def test_set_inlets_modes(self):
        pass

    def test_set_outlets_modes(self):
        pass

    def test_set_residual_modes(self):
        pass

    def test_set_boundary_conditions_clear(self):
        pass
