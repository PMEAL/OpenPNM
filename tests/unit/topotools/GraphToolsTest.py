import openpnm as op
import numpy as np
import pytest
from openpnm import topotools


class GraphToolsTest:

    def setup_class(self):
        r"""
        Make a nice and simple network to test stuff

        The Network:

        ::

            3 ― 4 ― 5
            | \ |
            0 ― 1 ― 2

            The Adjacency Matrix:    The Enumreated Adjacency Matrix:

              | 0  1  2  3  4  5       | 0  1  2  3  4  5
            ――|――――――――――――――――――    ――|――――――――――――――――――
            0 | -  1     1           0 | -  0     1
            1 | 1  -  1  1  1        1 | 0  -  3  2  5
            2 |    1  -              2 |    3  -
            3 | 1  1     -  1        3 | 1  2     -  4
            4 |    1     1  -  1     4 |    5     4  -  6
            5 |             1  -     5 |             6  -

            The Incidence Matrix:    The Enumerated Incidence Matrix

              | 0  1  2  3  4  5       | 0  1  2  3  4  5
            ――|――――――――――――――――――    ――|――――――――――――――――――
            0 | 1  1                 0 | 1  0
            1 | 1        1           1 | 3        0
            2 |    1     1           2 |    3     1
            3 |    1  1              3 |    2  1
            4 |          1  1        4 |          4  3
            5 |    1        1        5 |    4        1
            6 |             1  1     6 |             5  4

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
                          [1, 3],
                          [1, 2],
                          [3, 4],
                          [1, 4],
                          [4, 5]])
        self.net['pore.coords'] = coords
        self.net['throat.conns'] = conns

    def teardown_class(self):
        self.ws.clear()

    def test_find_connected_sites(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='coo')
        Ps = topotools.find_connected_sites(bonds=[0], am=am, flatten=True)
        assert np.all(Ps == [0, 1])
        Ps = topotools.find_connected_sites(bonds=[1], am=am, flatten=True)
        assert np.all(Ps == [0, 3])
        Ps = topotools.find_connected_sites(bonds=[2], am=am, flatten=True)
        assert np.all(Ps == [1, 3])
        Ps = topotools.find_connected_sites(bonds=[3], am=am, flatten=True)
        assert np.all(Ps == [1, 2])
        Ps = topotools.find_connected_sites(bonds=[4], am=am, flatten=True)
        assert np.all(Ps == [3, 4])
        Ps = topotools.find_connected_sites(bonds=[5], am=am, flatten=True)
        assert np.all(Ps == [1, 4])
        Ps = topotools.find_connected_sites(bonds=[6], am=am, flatten=True)
        assert np.all(Ps == [4, 5])

    def test_find_connected_sites_single(self):
        am = self.net.create_adjacency_matrix(fmt='coo')
        Ps = topotools.find_connected_sites(bonds=0, am=am, flatten=True)
        assert np.all(Ps == [0, 1])

    def test_find_connected_sites_fmt_not_coo(self):
        am = self.net.create_adjacency_matrix(fmt='csr')
        with pytest.raises(Exception):
            topotools.find_connected_sites(bonds=[0], am=am, flatten=True)

    def test_find_connected_sites_with_logic_flattened(self):
        am = self.net.get_adjacency_matrix(fmt='coo')
        a = topotools.find_connected_sites(bonds=[0, 2], am=am, logic='or',
                                           flatten=True)
        assert np.all(a == [0, 1, 3])
        a = topotools.find_connected_sites(bonds=[0, 2], am=am, logic='and',
                                           flatten=True)
        assert np.all(a == [])
        a = topotools.find_connected_sites(bonds=[0, 2], am=am, logic='xor',
                                           flatten=True)
        assert np.all(a == [0, 3])
        a = topotools.find_connected_sites(bonds=[0, 2], am=am, logic='xnor',
                                           flatten=True)
        assert np.all(a == [1])
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(bonds=[0], am=am, logic='nor')
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(bonds=[0], am=am, logic='nand')

    def test_find_connecting_bonds(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='dok')
        T = topotools.find_connecting_bonds(sites=[0, 1], am=am)
        assert np.all(T == [0])
        T = topotools.find_connecting_bonds(sites=[0, 3], am=am)
        assert np.all(T == [1])
        T = topotools.find_connecting_bonds(sites=[1, 3], am=am)
        assert np.all(T == [2])
        T = topotools.find_connecting_bonds(sites=[1, 2], am=am)
        assert np.all(T == [3])
        T = topotools.find_connecting_bonds(sites=[3, 4], am=am)
        assert np.all(T == [4])
        T = topotools.find_connecting_bonds(sites=[1, 4], am=am)
        assert np.all(T == [5])
        T = topotools.find_connecting_bonds(sites=[4, 5], am=am)
        assert np.all(T == [6])

    def test_find_connecting_bonds_fmt_not_dok(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='csr')
        T = topotools.find_connecting_bonds(sites=[0, 1], am=am)
        assert np.all(T == [0])
        T = topotools.find_connecting_bonds(sites=[0, 3], am=am)
        assert np.all(T == [1])
        T = topotools.find_connecting_bonds(sites=[1, 3], am=am)
        assert np.all(T == [2])
        T = topotools.find_connecting_bonds(sites=[1, 2], am=am)
        assert np.all(T == [3])
        T = topotools.find_connecting_bonds(sites=[3, 4], am=am)
        assert np.all(T == [4])
        T = topotools.find_connecting_bonds(sites=[1, 4], am=am)
        assert np.all(T == [5])
        T = topotools.find_connecting_bonds(sites=[4, 5], am=am)
        assert np.all(T == [6])

    def test_find_connecting_bonds_multiple_sites(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='dok')
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 3]], am=am)
        assert np.all(T == [0, 1])
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 3], [4, 5]],
                                            am=am)
        assert np.all(T == [0, 1, 6])

    def test_find_connecting_bonds_no_sites(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='dok')
        T = topotools.find_connecting_bonds(sites=[], am=am)
        assert T == []
        T = topotools.find_connecting_bonds(sites=[[], []], am=am)
        assert T == []

    def test_find_connecting_bonds_nonexistant_connections(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='dok')
        T = topotools.find_connecting_bonds(sites=[0, 5], am=am)
        assert np.all(T == [None])
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 5]], am=am)
        assert np.all(T == [0, None])

    def test_find_connecting_bonds_multiple_sites_fmt_not_dok(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='csr')
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 3]], am=am)
        assert np.all(T == [0, 1])
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 3], [4, 5]],
                                            am=am)
        assert np.all(T == [0, 1, 6])

    def test_find_neighbor_bonds(self):
        im = self.net.create_incidence_matrix(fmt='lil')
        Ts = topotools.find_neighbor_bonds(sites=[0], im=im)
        assert np.all(Ts == [0, 1])
        Ts = topotools.find_neighbor_bonds(sites=[1], im=im)
        assert np.all(Ts == [0, 2, 3, 5])
        Ts = topotools.find_neighbor_bonds(sites=[2], im=im)
        assert np.all(Ts == [3])
        Ts = topotools.find_neighbor_bonds(sites=[3], im=im)
        assert np.all(Ts == [1, 2, 4])
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
        assert np.all(Ts == [0, 2, 3, 5])
        Ts = topotools.find_neighbor_bonds(sites=[2], im=im)
        assert np.all(Ts == [3])
        Ts = topotools.find_neighbor_bonds(sites=[3], im=im)
        assert np.all(Ts == [1, 2, 4])
        Ts = topotools.find_neighbor_bonds(sites=[4], im=im)
        assert np.all(Ts == [4, 5, 6])
        Ts = topotools.find_neighbor_bonds(sites=[5], im=im)
        assert np.all(Ts == [6])

    def test_find_neighbor_bonds_with_logic(self):
        im = self.net.get_incidence_matrix(fmt='lil')
        a = topotools.find_neighbor_bonds(sites=[0, 2], im=im, logic='or')
        assert np.all(a == [0, 1, 3])
        a = topotools.find_neighbor_bonds(sites=[0, 2], im=im, logic='and')
        assert np.all(a == [])
        a = topotools.find_neighbor_bonds(sites=[0], im=im, logic='xor')
        assert np.all(a == [0, 1])
        a = topotools.find_neighbor_bonds(sites=[0], im=im, logic='xnor')
        assert np.all(a == [])
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(sites=[0], im=im, logic='nor')
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(sites=[0], im=im, logic='nand')

    def test_find_neighbor_bonds_with_am_and_logic(self):
        am = self.net.get_adjacency_matrix(fmt='coo')
        im = self.net.get_incidence_matrix(fmt='coo')
        Ts1 = topotools.find_neighbor_bonds(sites=[1], am=am,
                                            flatten=True, logic='or')
        Ts2 = topotools.find_neighbor_bonds(sites=[1], im=im,
                                            flatten=True, logic='or')
        assert np.all(Ts1 == Ts2)
        Ts1 = topotools.find_neighbor_bonds(sites=[1], am=am,
                                            flatten=True, logic='xor')
        Ts2 = topotools.find_neighbor_bonds(sites=[1], im=im,
                                            flatten=True, logic='xor')
        assert np.all(Ts1 == Ts2)
        Ts1 = topotools.find_neighbor_bonds(sites=[1], am=am,
                                            flatten=True, logic='xnor')
        Ts2 = topotools.find_neighbor_bonds(sites=[1], im=im,
                                            flatten=True, logic='xnor')
        assert np.all(Ts1 == Ts2)

    def test_find_neighbor_bonds_with_am_exceptions(self):
        am = self.net.get_adjacency_matrix(fmt='coo')
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(sites=[1], am=am, flatten=True,
                                          logic='intersection')
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(sites=[1], am=am, flatten=False,
                                          logic='or')

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

    def test_find_neighbor_sites_unflattened_or(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2], am=am, logic='or',
                                           flatten=False, include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [1, 3])
        assert np.all(Ps[1] == [0, 2, 3, 4])
        assert np.all(Ps[2] == [1])

    def test_find_neighbor_sites_unflattened_xor(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2], am=am,
                                           logic='xor',
                                           flatten=False, include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [])
        assert np.all(Ps[1] == [0, 2, 4])
        assert np.all(Ps[2] == [])

    def test_find_neighbor_sites_unflattened_and(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 2], am=am, logic='and',
                                           flatten=False, include_input=True)
        assert len(Ps) == 2
        assert np.all(Ps[0] == [1])
        assert np.all(Ps[1] == [1])

    def test_find_neighbor_sites_unflattened_xnor(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2], am=am,
                                           logic='xnor',
                                           flatten=False, include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [1, 3])
        assert np.all(Ps[1] == [3])
        assert np.all(Ps[2] == [1])

    def test_find_neighbor_sites_unflattened_xor_empty_set(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2], am=am, logic='or',
                                           flatten=False, include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [1, 3])
        assert np.all(Ps[1] == [0, 2, 3, 4])
        assert np.all(Ps[2] == [1])

    def test_find_neighbor_sites_unflattened_and_empty_set(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2], am=am, logic='and',
                                           flatten=False, include_input=False)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [])
        assert np.all(Ps[1] == [])
        assert np.all(Ps[2] == [])

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

    def test_find_neighbor_sites_flattened_with_logic(self):
        am = self.net.create_adjacency_matrix(fmt='lil')
        Ps = topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                           logic='or')
        assert np.all(Ps == [2, 3, 4])
        Ps = topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                           logic='and')
        assert np.all(Ps == [3])
        Ps = topotools.find_neighbor_sites(sites=[0, 1], am=am, flatten=True,
                                           logic='xor')
        assert np.all(Ps == [2, 4])
        Ps = topotools.find_neighbor_sites(sites=[0, 2], am=am, flatten=True,
                                           logic='xnor')
        assert np.all(Ps == [1])

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

    def test_find_complement(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='coo')
        a = [0, 1, 2]
        b = topotools.find_complement(sites=a, am=am)
        assert set(a).isdisjoint(b)
        a = [0, 1, 2]
        b = topotools.find_complement(bonds=a, am=am)
        assert set(a).isdisjoint(b)
        with pytest.raises(Exception):
            topotools.find_complement(am=am)
        with pytest.raises(Exception):
            topotools.find_complement(am=am, sites=a, bonds=a)

    def test_find_complement_asmask(self):
        am = self.net.create_adjacency_matrix(weights=self.net.Ts, fmt='coo')
        a = [0, 1, 2]
        b = topotools.find_complement(sites=a, am=am, asmask=True)
        assert len(b) == self.net.Np

if __name__ == '__main__':

    t = GraphToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
