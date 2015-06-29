import OpenPNM
import matplotlib as mpl


class PlotsTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        water = OpenPNM.Phases.Water(network=self.net)
        OpenPNM.Physics.Standard(network=self.net,
                                 phase=water,
                                 geometry=self.geo)

    def test_distributions(self):
        a = OpenPNM.Postprocessing.Plots.distributions(self.net)
        assert isinstance(a, mpl.figure.Figure)

    def test_pore_size_distribution(self):
        a = OpenPNM.Postprocessing.Plots.pore_size_distribution(self.net)
        assert isinstance(a, mpl.figure.Figure)

    def test_porosity_profile(self):
        a = OpenPNM.Postprocessing.Plots.porosity_profile(self.net, axis=0)
        assert isinstance(a, mpl.figure.Figure)
