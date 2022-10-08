import numpy as np
import pytest
from numpy.testing import assert_allclose
from openpnm._skgraph.generators import cubic
from openpnm._skgraph import queries
from openpnm._skgraph import tools


class SKGRQueriesTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_find_connected_nodes_undirected(self):
        g = cubic(shape=[3, 2, 1])
        c = queries.find_connected_nodes(network=g, inds=[0, 2, 4], logic='or')
        assert np.all(c == [0, 1, 3, 4, 5])
        c = queries.find_connected_nodes(network=g, inds=[0, 2, 4], logic='xor')
        assert np.all(c == [0, 3, 4, 5])
        c = queries.find_connected_nodes(network=g, inds=[0, 2, 4], logic='xnor')
        assert np.all(c == [1])
        c = queries.find_connected_nodes(network=g, inds=[0, 2, 4], logic='and')
        assert np.all(c == [])

    def test_find_connected_nodes_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        with pytest.raises(Exception):
            _ = queries.find_connected_nodes(network=g, inds=[0, 2, 4])

    def test_find_neighbor_edges_undirected(self):
        g = cubic(shape=[3, 2, 1])
        c = queries.find_neighbor_edges(network=g, inds=[0, 2, 4], logic='or')
        assert np.all(c == [0, 1, 2, 3, 5])
        c = queries.find_neighbor_edges(network=g, inds=[0, 2, 4], logic='xor')
        assert np.all(c == [0, 1, 2])
        c = queries.find_neighbor_edges(network=g, inds=[0, 2, 4], logic='xnor')
        assert np.all(c == [3, 5])
        with pytest.raises(Exception):
            _ = queries.find_neighbor_edges(network=g, inds=[0, 2, 4], logic='and')

    def test_find_neighbor_edges_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        c = queries.find_neighbor_edges(network=g, inds=[0, 1], logic='or')
        assert np.all(c == [0, 1, 3, 4])
        c = queries.find_neighbor_edges(network=g, inds=[0], logic='or')
        assert np.all(c == [0, 1, 3])
        c = queries.find_neighbor_edges(network=g, inds=[0, 1, 4], logic='xor')
        assert np.all(c == [2, 3, 4, 5])
        c = queries.find_neighbor_edges(network=g, inds=[0, 1, 4], logic='xnor')
        assert np.all(c == [0, 1])
        with pytest.raises(Exception):
            _ = queries.find_neighbor_edges(network=g, inds=[0, 2, 4], logic='and')

    def test_find_neighbor_nodes_undirected(self):
        g = cubic(shape=[3, 2, 1])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 2], logic='or')
        assert np.all(c == [1, 3, 4])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 2], logic='xor')
        assert np.all(c == [1, 3, 4])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 2], logic='xnor')
        assert np.all(c == [1])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 2], logic='and')
        assert np.all(c == [1])

    def test_find_neighbor_nodes_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [3, 2]
        c = queries.find_neighbor_nodes(network=g, inds=[0, 2], logic='or')
        assert np.all(c == [1, 4])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 2], logic='xor')
        assert np.all(c == [1, 4])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 3], logic='xnor')
        assert np.all(c == [2])
        c = queries.find_neighbor_nodes(network=g, inds=[0, 3], logic='and')
        assert np.all(c == [2])

    def test_find_connecting_edges_undirected(self):
        g = cubic(shape=[3, 2, 1])
        c = queries.find_connecting_edges(inds=[0, 1], network=g)
        assert np.all(c == [0])
        c = queries.find_connecting_edges(inds=[[0, 1]], network=g)
        assert np.all(c == [0])
        c = queries.find_connecting_edges(inds=[[0, 1], [0, 2]], network=g)
        assert np.all(c == [0, 3])
        c = queries.find_connecting_edges(inds=[[0, 1], [0, 2], [0, 5]], network=g)
        assert np.all(c[:2] == [0., 3.])
        assert np.isnan(c[2])

    def test_find_connecting_edges_directed(self):
        g = cubic(shape=[3, 2, 1])
        edge_prefix = tools.get_edge_prefix(g)
        g[edge_prefix+'.conns'][1, :] = [1, 0]
        c = queries.find_connecting_edges(inds=[0, 1], network=g)
        assert np.all(c == [0])
        c = queries.find_connecting_edges(inds=[1, 0], network=g)
        assert np.all(c == [1])
        c = queries.find_connecting_edges(inds=[[0, 1], [1, 0]], network=g)
        assert np.all(c == [0, 1])
        c = queries.find_connecting_edges(inds=[[0, 1], [1, 0], [1, 5]], network=g)
        assert np.all(c[:2] == [0., 1.])
        assert np.isnan(c[2])

    def test_find_coordination_undirected(self):
        g = cubic(shape=[3, 2, 1])
        z = queries.find_coordination(g)
        assert np.all(z == [2, 2, 3, 3, 2, 2])
        z = queries.find_coordination(g, nodes=[0, 2])
        assert np.all(z == [2, 3])

    def test_find_coordination_directed_symmetrically(self):
        g = cubic(shape=[3, 2, 1])
        edge_prefix = tools.get_edge_prefix(g)
        conns = g[edge_prefix+'.conns']
        g[edge_prefix+'.conns'] = np.vstack((conns, np.fliplr(conns)))
        z = queries.find_coordination(g)
        assert np.all(z == [2, 2, 3, 3, 2, 2])
        z = queries.find_coordination(g, nodes=[0, 2])
        assert np.all(z == [2, 3])

    def test_find_coordination_directed_nonsymmetrically(self):
        g = cubic(shape=[3, 2, 1])
        edge_prefix = tools.get_edge_prefix(g)
        g[edge_prefix+'.conns'][1, :] = [5, 1]
        z = queries.find_coordination(g)
        assert np.all(z == [2, 1, 1, 1, 1, 1])
        z = queries.find_coordination(g, nodes=[0, 2])
        assert np.all(z == [2, 1])

    def test_filter_by_z_undirected(self):
        g = cubic(shape=[3, 2, 1])
        z = queries.filter_by_z(network=g, inds=[0, 2, 3], z=2)
        assert np.all(z == [0])
        z = queries.filter_by_z(network=g, inds=[0, 2, 3], z=3)
        assert np.all(z == [2, 3])
        z = queries.filter_by_z(network=g, inds=[0, 2, 3], z=4)
        assert np.all(z == [])

    def test_filter_by_z_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        z = queries.filter_by_z(network=g, inds=[0, 2, 3], z=3)
        assert z.size == 0
        z = queries.filter_by_z(network=g, inds=[0, 2, 3], z=2)
        assert np.all(z == [0])

    def test_find_common_edges_undirected(self):
        g = cubic(shape=[3, 2, 1])
        c = queries.find_common_edges(network=g, inds_1=[0, 1], inds_2=[4, 5])
        assert len(c) == 0
        c = queries.find_common_edges(network=g, inds_1=[0, 1], inds_2=[3, 4])
        assert c == np.array([4])
        c = queries.find_common_edges(network=g, inds_1=[1, 2], inds_2=[3, 4])
        assert np.all(c == np.array([1, 4, 5]))

    def test_find_common_edges_intersecting_inds(self):
        g = cubic(shape=[3, 2, 1])
        with pytest.raises(Exception):
            _ = queries.find_common_edges(network=g, inds_1=[0, 1], inds_2=[0, 1])

    def test_find_common_edges_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        with pytest.raises(Exception):
            _ = queries.find_common_edges(network=g, inds_1=[0, 1], inds_2=[4, 5])

    def test_find_complementary_edges_undirected(self):
        g = cubic(shape=[2, 2, 1])
        c = queries.find_complementary_edges(inds=[0, 2], network=g)
        assert np.all(c == np.array([1, 3]))
        c = queries.find_complementary_edges(inds=[0, 2], network=g, asmask=True)
        assert np.all(c == np.array([False, True, False, True]))

    def test_find_complementary_edges_directed(self):
        g = cubic(shape=[3, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        c = queries.find_complementary_edges(inds=[0, 1, 2], network=g)
        assert np.all(c == np.array([3, 4, 5, 6]))
        c = queries.find_complementary_edges(inds=[0, 1, 2], network=g, asmask=True)
        assert np.all(c == np.array([False, False, False, True, True, True, True]))

    def test_find_complementary_nodes(self):
        g = cubic(shape=[2, 2, 1])
        c = queries.find_complementary_nodes(inds=[0, 2], network=g)
        assert np.all(c == np.array([1, 3]))
        c = queries.find_complementary_nodes(inds=[0, 2], network=g, asmask=True)
        assert np.all(c == np.array([False, True, False, True]))

    def test_find_complementary_nodes_directed(self):
        # Network directedness should have not impact on this, so using same
        # test as above
        g = cubic(shape=[2, 2, 1])
        g['edge.conns'][1, :] = [1, 0]
        c = queries.find_complementary_nodes(inds=[0, 2], network=g)
        assert np.all(c == np.array([1, 3]))
        c = queries.find_complementary_nodes(inds=[0, 2], network=g, asmask=True)
        assert np.all(c == np.array([False, True, False, True]))

    def test_find_path_undirected(self):
        g = cubic(shape=[4, 4, 1])
        p = queries.find_path(network=g, pairs=[[0, 11], [11, 0]])
        a = np.unique(g['edge.conns'][p['edge_paths'][0]])
        b = np.sort(p['node_paths'][0])
        assert np.all(a == b)
        a = np.unique(g['edge.conns'][p['edge_paths'][1]])
        b = np.sort(p['node_paths'][1])
        assert np.all(a == b)

    def test_find_path_undirected_neighboring_pores(self):
        g = cubic(shape=[4, 4, 1])
        p = queries.find_path(network=g, pairs=[[0, 1], [1, 0]])
        assert np.all(p['node_paths'][0] == [0, 1])
        assert np.all(p['node_paths'][1] == [1, 0])
        assert p['edge_paths'][0] == [0]
        assert p['edge_paths'][1] == [0]

    def test_find_path_directed(self):
        g = cubic(shape=[4, 4, 1])
        g['edge.conns'][0, :] = [1, 0]
        p = queries.find_path(network=g, pairs=[[0, 1], [1, 0]])
        assert np.all(p['node_paths'][0] == [])
        assert np.all(p['node_paths'][1] == [1, 0])
        assert p['edge_paths'][0] == []
        assert p['edge_paths'][1] == [0]


if __name__ == '__main__':
    t = SKGRQueriesTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
