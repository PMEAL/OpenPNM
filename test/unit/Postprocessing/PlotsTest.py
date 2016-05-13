import OpenPNM
import matplotlib as mpl


class PlotsTest:
    def setup_class(self):
        mgr = OpenPNM.Base.Workspace()
        mgr.clear()
        self.net = OpenPNM.Network.Cubic(shape=[10, 10, 10])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.water = OpenPNM.Phases.Water(network=self.net)
        self.phys = OpenPNM.Physics.Standard(network=self.net,
                                             phase=self.water,
                                             geometry=self.geo)
        self.IP = OpenPNM.Algorithms.InvasionPercolation(network=self.net)
        self.IP.setup(phase=self.water)
        self.IP.set_inlets(self.net.pores('top'))
        self.IP.run()

    def test_distributions(self):
        a = OpenPNM.Postprocessing.Plots.distributions(self.net)
        assert isinstance(a, mpl.figure.Figure)

    def test_pore_size_distribution(self):
        a = OpenPNM.Postprocessing.Plots.pore_size_distribution(self.net)
        assert isinstance(a, mpl.figure.Figure)

    def test_porosity_profile_axis_0(self):
        a = OpenPNM.Postprocessing.Plots.porosity_profile(self.net, axis=0)
        assert isinstance(a, mpl.figure.Figure)

    def test_porosity_profile_axis_2(self):
        a = OpenPNM.Postprocessing.Plots.porosity_profile(self.net, axis=1)
        assert isinstance(a, mpl.figure.Figure)

    def test_porosity_profile_axis_3(self):
        a = OpenPNM.Postprocessing.Plots.porosity_profile(self.net, axis=2)
        assert isinstance(a, mpl.figure.Figure)

#    def test_profiles(self):
#        vals = self.geo['pore.diameter']
#        a = OpenPNM.Postprocessing.Plots.profiles(self.net, values=vals)
#        assert isinstance(a, mpl.figure.Figure)
#
#    def test_saturation_profile(self):
#        self.water['pore.occupancy'] = self.IP['pore.invaded'] < 500
#        self.water['throat.occupancy'] = self.IP['throat.invaded'] < 500
#        a = OpenPNM.Postprocessing.Plots.saturation_profile(self.net,
#                                                            phase=self.water)
#        assert isinstance(a, mpl.figure.Figure)
