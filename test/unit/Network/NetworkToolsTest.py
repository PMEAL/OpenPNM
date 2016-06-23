import OpenPNM as op
import scipy as sp
import matplotlib as mpl


class NetworkToolsTest:
    def setup_class(self):
        self.net = op.Network.Cubic(shape=[5, 5, 5], spacing=1)
        self.net['pore.diameter'] = sp.rand(self.net.Np)

    def test_plot_network(self):
        a = op.Network.tools.plot_network(self.net)
        assert isinstance(a, mpl.figure.Figure)
