import OpenPNM as op
import scipy as sp
import matplotlib as mpl


class NetworkToolsTest:
    def setup_class(self):
        self.net = op.Network.Cubic(shape=[5, 5, 5], spacing=1)
        self.net['pore.diameter'] = sp.rand(self.net.Np)

    def test_plot_topology(self):
        a = op.Network.tools.plot_topology(self.net)
        assert isinstance(a, mpl.figure.Figure)
        b = op.Network.tools.plot_topology(network=self.net, fig=a, c='r')
        assert b is a
        c = op.Network.tools.plot_topology(network=self.net, fig=b,
                                           throats=[1, 2, 3], c='b')
        assert c is b

    def test_plot_coordinates(self):
        a = op.Network.tools.plot_coordinates(self.net)
        assert isinstance(a, mpl.figure.Figure)
        b = op.Network.tools.plot_coordinates(network=self.net, fig=a, c='r')
        assert b is a
        c = op.Network.tools.plot_coordinates(network=self.net, fig=b,
                                              pores=[1, 2, 3], c='b', s=50)
        assert c is b
