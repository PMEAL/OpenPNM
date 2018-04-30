import openpnm as op
import numpy as np
import pytest
from openpnm import topotools


class GraphToolsTest:

    def setup_class(self):
        r"""
        Make a nice and simple network to test stuff

        3 ― 4 ― 5
        | \ |
        0 ― 1 ― 2

        """
        self.ws = op.Workspace()
        self.net = op.network.GenericNetwork(Np=6, Nt=7)
        coords = np.array([[0, 0, 0],
                           [1, 0, 0],
                           [2, 0, 0],
                           [0, 1, 0],
                           [1, 1, 0],
                           [2, 1, 0]])
        conns = np.array([[0, 1],
                          [0, 3],
                          [1, 2],
                          [1, 3],
                          [1, 4],
                          [3, 4],
                          [4, 5]])
        self.net['pore.coords'] = coords
        self.net['throat.conns'] = conns

    def teardown_class(self):
        self.ws.clear()

    def test_find_neighbor_bonds(self):
        im = self.net.create_incidence_matrix(fmt='lil')
        Ts = topotools.find_neighbor_bonds(sites=[0], im=im)
        assert np.all(Ts == [0, 1])
        Ts = topotools.find_neighbor_bonds(sites=[1], im=im)
        assert np.all(Ts == [0, 2, 3, 4])
        Ts = topotools.find_neighbor_bonds(sites=[2], im=im)
        assert np.all(Ts == [2])
        Ts = topotools.find_neighbor_bonds(sites=[3], im=im)
        assert np.all(Ts == [1, 3, 5])
        Ts = topotools.find_neighbor_bonds(sites=[4], im=im)
        assert np.all(Ts == [4, 5, 6])
        Ts = topotools.find_neighbor_bonds(sites=[5], im=im)
        assert np.all(Ts == [6])

    def test_find_neighbor_bonds_single(self):
        im = self.net.create_incidence_matrix(fmt='lil')
        Ts = topotools.find_neighbor_bonds(sites=0, im=im)
        assert np.all(Ts == [0, 1])

    def test_find_neighbor_bonds_fmt_not_lil(self):
        # Make sure results are correct even if converting to lil internally
        im = self.net.create_incidence_matrix(fmt='coo')
        Ts = topotools.find_neighbor_bonds(sites=[0], im=im)
        assert np.all(Ts == [0, 1])
        Ts = topotools.find_neighbor_bonds(sites=[1], im=im)
        assert np.all(Ts == [0, 2, 3, 4])
        Ts = topotools.find_neighbor_bonds(sites=[2], im=im)
        assert np.all(Ts == [2])
        Ts = topotools.find_neighbor_bonds(sites=[3], im=im)
        assert np.all(Ts == [1, 3, 5])
        Ts = topotools.find_neighbor_bonds(sites=[4], im=im)
        assert np.all(Ts == [4, 5, 6])
        Ts = topotools.find_neighbor_bonds(sites=[5], im=im)
        assert np.all(Ts == [6])

    def test_find_neighbor_sites(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0], am=am)
        assert np.all(Ps == [1, 3])
        Ps = topotools.find_neighbor_sites(sites=[1], am=am)
        assert np.all(Ps == [0, 2, 3, 4])
        Ps = topotools.find_neighbor_sites(sites=[2], am=am)
        assert np.all(Ps == [1])
        Ps = topotools.find_neighbor_sites(sites=[3], am=am)
        assert np.all(Ps == [0, 1, 4])
        Ps = topotools.find_neighbor_sites(sites=[4], am=am)
        assert np.all(Ps == [1, 3, 5])
        Ps = topotools.find_neighbor_sites(sites=[5], am=am)
        assert np.all(Ps == [4])

    def test_find_neighbor_sites_single(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=0, am=am)
        assert np.all(Ps == [1, 3])

    def test_find_neighbor_sites_fmt_not_lil(self):
        # Make sure results are correct even if converting to lil internally
        am = self.net.create_adjacency_matrix(fmt='csr')
        Ps = topotools.find_neighbor_sites(sites=[0], am=am)
        assert np.all(Ps == [1, 3])
        Ps = topotools.find_neighbor_sites(sites=[1], am=am)
        assert np.all(Ps == [0, 2, 3, 4])
        Ps = topotools.find_neighbor_sites(sites=[2], am=am)
        assert np.all(Ps == [1])
        Ps = topotools.find_neighbor_sites(sites=[3], am=am)
        assert np.all(Ps == [0, 1, 4])
        Ps = topotools.find_neighbor_sites(sites=[4], am=am)
        assert np.all(Ps == [1, 3, 5])
        Ps = topotools.find_neighbor_sites(sites=[5], am=am)
        assert np.all(Ps == [4])

    def test_find_neighbor_sites_with_logic(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                           logic='union')
        assert np.all(Ps == [2, 3, 4])
        Ps = topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                           logic='intersection')
        assert np.all(Ps == [3])
        Ps = topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                           logic='exclusive_or')
        assert np.all(Ps == [2, 4])

    def test_find_neighbor_sites_with_bad_logic(self):
        with pytest.raises(Exception):
            am = self.net.create_adjacency_matrix(fmt='lil')
            topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                          logic='foobar')

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

    def test_find_clusters_sites(self):
        net = op.network.Cubic(shape=[10, 10, 10])
        net['pore.seed'] = np.random.rand(net.Np)
        net['throat.seed'] = np.random.rand(net.Nt)
        clusters = topotools.find_clusters(network=net,
                                           mask=net['pore.seed'] < 0.5)
        assert len(clusters[0]) == net.Np
        assert len(clusters[1]) == net.Nt

    def test_find_clusters_bonds(self):
        net = op.network.Cubic(shape=[10, 10, 10])
        net['pore.seed'] = np.random.rand(net.Np)
        net['throat.seed'] = np.random.rand(net.Nt)
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
