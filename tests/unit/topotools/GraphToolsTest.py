import openpnm as op
import numpy as np
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
        # Now test non-coo AM's
        am = net.create_adjacency_matrix(triu=False, fmt='lil')
        assert not topotools.istriu(am)
        am = net.create_adjacency_matrix(triu=True, fmt='csr')
        assert topotools.istriu(am)
        am = am.T
        assert not topotools.istriu(am)
        # Now test non-triangular AM
        im = net.create_incidence_matrix()
        assert not topotools.istriu(im)

    def test_istril(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istril(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.istril(am)
        am = am.T
        assert topotools.istril(am)
        # Now test non-coo AM's
        am = net.create_adjacency_matrix(triu=False, fmt='lil')
        assert not topotools.istril(am)
        am = net.create_adjacency_matrix(triu=True, fmt='csr')
        assert not topotools.istril(am)
        am = am.T
        assert topotools.istril(am)
        # Now test non-triangular AM
        im = net.create_incidence_matrix()
        assert not topotools.istril(im)

    def test_istriangular(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istriangular(am)
        am = net.create_adjacency_matrix(triu=True)
        assert topotools.istriangular(am)
        am = am.T
        assert topotools.istriangular(am)
        # Now test non-coo AM's
        am = net.create_adjacency_matrix(triu=False, fmt='lil')
        assert not topotools.istriangular(am)
        am = net.create_adjacency_matrix(triu=True, fmt='csr')
        assert topotools.istriangular(am)
        am = am.T
        assert topotools.istriangular(am)
        # Now test non-triangular AM
        im = net.create_incidence_matrix()
        assert not topotools.istriangular(im)

    def test_issymmetric(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert topotools.issymmetric(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.issymmetric(am)
        am = am.T
        assert not topotools.issymmetric(am)
        # Now test non-coo AM's
        am = net.create_adjacency_matrix(triu=False)
        assert topotools.issymmetric(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.issymmetric(am)
        am = am.T
        assert not topotools.issymmetric(am)
        # Now test non-triangular AM
        im = net.create_incidence_matrix()
        assert not topotools.issymmetric(im)

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
