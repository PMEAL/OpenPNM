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

            The Adjacency Matrix:    The Enumerated Adjacency Matrix:

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
        self.net = op.network.Network()
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
        Ps = topotools.find_connected_sites(bonds=[0],network=self.net, flatten=True)
        assert np.all(Ps == [0, 1])
        Ps = topotools.find_connected_sites(bonds=[1], network=self.net, flatten=True)
        assert np.all(Ps == [0, 3])
        Ps = topotools.find_connected_sites(bonds=[2], network=self.net, flatten=True)
        assert np.all(Ps == [1, 3])
        Ps = topotools.find_connected_sites(bonds=[3], network=self.net, flatten=True)
        assert np.all(Ps == [1, 2])
        Ps = topotools.find_connected_sites(bonds=[4], network=self.net, flatten=True)
        assert np.all(Ps == [3, 4])
        Ps = topotools.find_connected_sites(bonds=[5], network=self.net, flatten=True)
        assert np.all(Ps == [1, 4])
        Ps = topotools.find_connected_sites(bonds=[6], network=self.net, flatten=True)
        assert np.all(Ps == [4, 5])

    def test_find_connected_sites_unflattened(self):
        Ps = topotools.find_connected_sites(bonds=[0, 5],
                                            network=self.net,
                                            flatten=False,
                                            logic="xnor")
        x = np.array([[np.nan, 1], [1, np.nan]])
        assert ((x == Ps) | (np.isnan(x) & np.isnan(Ps))).all()
        Ps = topotools.find_connected_sites(bonds=[0, 1, 5],
                                            network=self.net,
                                            logic="and",
                                            flatten=False)
        assert np.allclose(Ps, [np.array([], dtype="int64")] * 3)

    def test_find_connected_sites_single(self):
        Ps = topotools.find_connected_sites(bonds=0, network=self.net, flatten=True)
        assert np.all(Ps == [0, 1])

    def test_find_connected_sites_with_logic_flattened(self):
        a = topotools.find_connected_sites(bonds=[0, 2],
                                           network=self.net,
                                           logic='or',
                                           flatten=True)
        assert np.all(a == [0, 1, 3])
        a = topotools.find_connected_sites(bonds=[0, 2],
                                           network=self.net,
                                           logic='and',
                                           flatten=True)
        assert np.all(a == [])
        a = topotools.find_connected_sites(bonds=[0, 2],
                                           network=self.net,
                                           logic='xor',
                                           flatten=True)
        assert np.all(a == [0, 3])
        a = topotools.find_connected_sites(bonds=[0, 2],
                                           network=self.net,
                                           logic='xnor',
                                           flatten=True)
        assert np.all(a == [1])
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(bonds=[0],
                                          network=self.net,
                                          logic='nor')
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(bonds=[0],
                                          network=self.net,
                                          logic='nand')

    def test_find_connecting_bonds(self):
        T = topotools.find_connecting_bonds(sites=[0, 1], network=self.net)
        assert np.all(T == [0])
        T = topotools.find_connecting_bonds(sites=[0, 3], network=self.net)
        assert np.all(T == [1])
        T = topotools.find_connecting_bonds(sites=[1, 3], network=self.net)
        assert np.all(T == [2])
        T = topotools.find_connecting_bonds(sites=[1, 2], network=self.net)
        assert np.all(T == [3])
        T = topotools.find_connecting_bonds(sites=[3, 4], network=self.net)
        assert np.all(T == [4])
        T = topotools.find_connecting_bonds(sites=[1, 4], network=self.net)
        assert np.all(T == [5])
        T = topotools.find_connecting_bonds(sites=[4, 5], network=self.net)
        assert np.all(T == [6])

    def test_find_connecting_bonds_multiple_sites(self):
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 3]], network=self.net)
        assert np.all(T == [0, 1])
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 3], [4, 5]],
                                            network=self.net)
        assert np.all(T == [0, 1, 6])

    def test_find_connecting_bonds_no_sites(self):
        T = topotools.find_connecting_bonds(sites=[], network=self.net)
        assert T == []
        T = topotools.find_connecting_bonds(sites=[[], []], network=self.net)
        assert T == []

    def test_find_connecting_bonds_nonexistant_connections(self):
        T = topotools.find_connecting_bonds(sites=[0, 5], network=self.net)
        assert np.isnan(T).sum() == 1
        T = topotools.find_connecting_bonds(sites=[[0, 1], [0, 5]], network=self.net)
        assert np.isnan(T).sum() == 1

    def test_find_neighbor_bonds(self):
        Ts = topotools.find_neighbor_bonds(sites=[0], network=self.net)
        assert np.all(Ts == [0, 1])
        Ts = topotools.find_neighbor_bonds(sites=[1], network=self.net)
        assert np.all(Ts == [0, 2, 3, 5])
        Ts = topotools.find_neighbor_bonds(sites=[2], network=self.net)
        assert np.all(Ts == [3])
        Ts = topotools.find_neighbor_bonds(sites=[3], network=self.net)
        assert np.all(Ts == [1, 2, 4])
        Ts = topotools.find_neighbor_bonds(sites=[4], network=self.net)
        assert np.all(Ts == [4, 5, 6])
        Ts = topotools.find_neighbor_bonds(sites=[5], network=self.net)
        assert np.all(Ts == [6])
        Ts = topotools.find_neighbor_bonds(sites=[], network=self.net)
        assert np.all(Ts == [])

    def test_find_neighbr_bonds_unflattened(self):
        Ts = topotools.find_neighbor_bonds(sites=[0, 1, 5],
                                           logic="and",
                                           network=self.net,
                                           flatten=False)
        assert np.allclose(Ts, [np.array([], dtype="int64")] * 3)

    def test_find_neighbor_bonds_given_am(self):
        Ts = topotools.find_neighbor_bonds(sites=[0], network=self.net)
        assert np.all(Ts == [0, 1])
        with pytest.raises(Exception):
            _ = topotools.find_neighbor_bonds(sites=[0], network=self.net,
                                              logic="unsupported_logic")

    def test_find_neighbor_bonds_missing_both_am_and_im(self):
        with pytest.raises(Exception):
            _ = topotools.find_neighbor_bonds(sites=[0])

    def test_find_neighbor_bonds_single(self):
        Ts = topotools.find_neighbor_bonds(sites=0, network=self.net)
        assert np.all(Ts == [0, 1])

    def test_find_neighbor_bonds_with_logic(self):
        a = topotools.find_neighbor_bonds(sites=[0, 2], network=self.net, logic='or')
        assert np.all(a == [0, 1, 3])
        a = topotools.find_neighbor_bonds(sites=[0], network=self.net, logic='xor')
        assert np.all(a == [0, 1])
        a = topotools.find_neighbor_bonds(sites=[0], network=self.net, logic='xnor')
        assert np.all(a == [])
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(sites=[0], network=self.net, logic='nor')
        with pytest.raises(Exception):
            topotools.find_neighbor_bonds(sites=[0], network=self.net, logic='nand')

    def test_find_neighbor_sites(self):
        Ps = topotools.find_neighbor_sites(sites=[0], network=self.net)
        assert np.all(Ps == [1, 3])
        Ps = topotools.find_neighbor_sites(sites=[1], network=self.net)
        assert np.all(Ps == [0, 2, 3, 4])
        Ps = topotools.find_neighbor_sites(sites=[2], network=self.net)
        assert np.all(Ps == [1])
        Ps = topotools.find_neighbor_sites(sites=[3], network=self.net)
        assert np.all(Ps == [0, 1, 4])
        Ps = topotools.find_neighbor_sites(sites=[4], network=self.net)
        assert np.all(Ps == [1, 3, 5])
        Ps = topotools.find_neighbor_sites(sites=[5], network=self.net)
        assert np.all(Ps == [4])
        Ps = topotools.find_neighbor_sites(sites=[], network=self.net)
        assert Ps == []

    def test_find_neighbor_sites_unsupported_logic(self):
        with pytest.raises(Exception):
            _ = topotools.find_neighbor_sites(sites=[0], network=self.net,
                                              logic="unsupported_logic")

    def test_find_neighbor_sites_single(self):
        Ps = topotools.find_neighbor_sites(sites=0, network=self.net)
        assert np.all(Ps == [1, 3])

    def test_find_neighbor_sites_unflattened_or(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2],
                                           network=self.net,
                                           logic='or',
                                           flatten=False,
                                           include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [1, 3])
        assert np.all(Ps[1] == [0, 2, 3, 4])
        assert np.all(Ps[2] == [1])

    def test_find_neighbor_sites_unflattened_xor(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2],
                                           network=self.net,
                                           logic='xor',
                                           flatten=False,
                                           include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [])
        assert np.all(Ps[1] == [0, 2, 4])
        assert np.all(Ps[2] == [])

    def test_find_neighbor_sites_unflattened_and(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 2],
                                           network=self.net,
                                           logic='and',
                                           flatten=False,
                                           include_input=True)
        assert len(Ps) == 2
        assert np.all(Ps[0] == [1])
        assert np.all(Ps[1] == [1])

    def test_find_neighbor_sites_unflattened_xnor(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2],
                                           network=self.net,
                                           logic='xnor',
                                           flatten=False,
                                           include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [1, 3])
        assert np.all(Ps[1] == [3])
        assert np.all(Ps[2] == [1])

    def test_find_neighbor_sites_unflattened_xor_empty_set(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2],
                                           network=self.net,
                                           logic='or',
                                           flatten=False,
                                           include_input=True)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [1, 3])
        assert np.all(Ps[1] == [0, 2, 3, 4])
        assert np.all(Ps[2] == [1])

    def test_find_neighbor_sites_unflattened_and_empty_set(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1, 2],
                                           network=self.net,
                                           logic='and',
                                           flatten=False,
                                           include_input=False)
        assert len(Ps) == 3
        assert np.all(Ps[0] == [])
        assert np.all(Ps[1] == [])
        assert np.all(Ps[2] == [])

    def test_find_neighbor_sites_flattened_with_logic(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1],
                                           network=self.net,
                                           flatten=True,
                                           logic='or')
        assert np.all(Ps == [2, 3, 4])
        Ps = topotools.find_neighbor_sites(sites=[0, 1],
                                           network=self.net,
                                           flatten=True,
                                           logic='and')
        assert np.all(Ps == [3])
        Ps = topotools.find_neighbor_sites(sites=[0, 1],
                                           network=self.net,
                                           flatten=True,
                                           logic='xor')
        assert np.all(Ps == [2, 4])
        Ps = topotools.find_neighbor_sites(sites=[0, 2],
                                           network=self.net,
                                           flatten=True,
                                           logic='xnor')
        assert np.all(Ps == [1])

    def test_find_neighbor_sites_with_bad_logic(self):
        with pytest.raises(Exception):
            topotools.find_neighbor_sites(sites=[0, 1],
                                          network=self.net,
                                          flatten=True,
                                          logic='foobar')

    def test_find_neighbor_sites_include_inputs(self):
        Ps = topotools.find_neighbor_sites(sites=[0, 1],
                                           network=self.net,
                                           flatten=True,
                                           logic='or',
                                           include_input=True)
        assert (Ps == [0, 1, 2, 3, 4]).all()
        Ps = topotools.find_neighbor_sites(sites=[0, 1],
                                           network=self.net,
                                           flatten=True,
                                           logic='or',
                                           include_input=False)
        assert (Ps == [2, 3, 4]).all()

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
        am = net.create_adjacency_matrix(triu=False, fmt="lil")
        assert topotools.issymmetric(am)
        am = net.create_adjacency_matrix(triu=True, fmt="lil")
        assert not topotools.issymmetric(am)
        am = am.T
        assert not topotools.issymmetric(am)
        # Now test non-triangular AM
        im = net.create_incidence_matrix()
        assert not topotools.issymmetric(im)

    def test_drop_sites(self):
        am = self.net.create_adjacency_matrix(fmt='coo')
        assert np.all(am.shape == (6, 6))
        assert am.col.max() == 5
        am, Ts = topotools.drop_sites(am, sites=[0])
        assert np.all(am.shape == (5, 5))
        assert am.col.max() == 4


if __name__ == '__main__':
    t = GraphToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
