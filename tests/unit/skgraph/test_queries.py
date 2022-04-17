import numpy as np
from numpy.testing import assert_allclose
from openpnm._skgraph.generators import cubic
from openpnm._skgraph import queries
from openpnm._skgraph import tools


class SKGRQueriesTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_find_coordination_undirected(self):
        g = cubic(shape=[3, 2, 1])
        z = queries.find_coordination(g)
        assert np.all(z == [2, 2, 3, 3, 2, 2])
        z = queries.find_coordination(g, nodes=[0, 2])
        assert np.all(z == [2, 3])

    def test_find_coordination_directed(self):
        g = cubic(shape=[3, 2, 1])
        edge_prefix = tools.get_edge_prefix(g)
        conns = g[edge_prefix+'.conns']
        g[edge_prefix+'.conns'] = np.vstack((conns, np.fliplr(conns)))
        z = queries.find_coordination(g)
        assert np.all(z == [2, 2, 3, 3, 2, 2])
        z = queries.find_coordination(g, nodes=[0, 2])
        assert np.all(z == [2, 3])

    def test_filter_by_z(self):
        g = cubic(shape=[3, 2, 1])
        z = queries.filter_by_z(g=g, inds=[0, 2, 3], z=2)
        assert np.all(z == [0])
        z = queries.filter_by_z(g=g, inds=[0, 2, 3], z=3)
        assert np.all(z == [2, 3])
        z = queries.filter_by_z(g=g, inds=[0, 2, 3], z=4)
        assert np.all(z == [])

    def test_find_common_edges(self):
        g = cubic(shape=[3, 2, 1])
        c = queries.find_common_edges(g=g, inds_1=[0, 1], inds_2=[4, 5])


if __name__ == '__main__':
    t = SKGRQueriesTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
