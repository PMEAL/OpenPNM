import openpnm as op
import numpy as np
from numpy.testing import assert_approx_equal
from openpnm import topotools


class GraphToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_istriu(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istriu(am)
        am = net.create_adjacency_matrix(triu=True)
        assert topotools.istriu(am)
        am = am.T
        assert not topotools.istriu(am)

    def test_istril(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istril(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.istril(am)
        am = am.T
        assert topotools.istril(am)

    def test_istriangular(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istriangular(am)
        am = net.create_adjacency_matrix(triu=True)
        assert topotools.istriangular(am)
        am = am.T
        assert topotools.istriangular(am)

    def test_issymmetric(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert topotools.issymmetric(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.issymmetric(am)
        am = am.T
        assert not topotools.issymmetric(am)

    def test_find_clusters(self):
        net = op.network.Cubic(shape=[10, 10, 10])
        net['pore.seed'] = np.random.rand(net.Np)
        net['throat.seed'] = np.random.rand(net.Nt)
        clusters = topotools.find_clusters(network=net,
                                           mask=net['pore.seed'] < 0.5)
        assert len(clusters[0]) == net.Np
        assert len(clusters[1]) == net.Nt
        clusters = topotools.find_clusters(network=net,
                                           mask=net['throat.seed'] < 0.5)
        assert len(clusters[0]) == net.Np
        assert len(clusters[1]) == net.Nt


if __name__ == '__main__':

    t = GraphToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
