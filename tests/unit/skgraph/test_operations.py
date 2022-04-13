import numpy as np
from openpnm._skgraph import generators as gen
from openpnm._skgraph import settings
from openpnm._skgraph.visualization import plot_edges, plot_nodes
import openpnm._skgraph.operations as ops
settings.node_prefix = 'node'
settings.edge_prefix = 'edge'


class SKGROperationsTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_trim_edges(self):
        net = gen.cubic([3, 3, 3])
        net['edge.label'] = np.ones(54, dtype=bool)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 54
        assert np.all(net['edge.conns'][0] == [0, 1])
        net = ops.trim_edges(net, 0)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 53
        assert net['edge.label'].shape[0] == 53
        assert np.all(net['edge.conns'][0] == [1, 2])
        net = ops.trim_edges(net, [20, 30])
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 51
        assert net['edge.label'].shape[0] == 51
        # ax = plot_edges(net['edge.conns'], net['node.coords'])
        # ax = plot_nodes(net['node.coords'], ax=ax)

    def test_trim_nodes(self):
        net = gen.cubic([4, 4, 1])
        net['edge.label'] = np.ones(24, dtype=bool)
        assert net['node.coords'].shape[0] == 16
        assert net['edge.conns'].shape[0] == 24
        assert np.all(net['edge.conns'][0] == [0, 1])
        net = ops.trim_nodes(net, 0)
        assert net['node.coords'].shape[0] == 15
        assert net['edge.conns'].shape[0] == 22
        assert net['edge.label'].shape[0] == 22
        assert np.all(net['edge.conns'][0] == [0, 1])
        net = ops.trim_nodes(net, [5, 10])
        assert net['node.coords'].shape[0] == 13
        assert net['edge.conns'].shape[0] == 15
        assert net['edge.label'].shape[0] == 15
        # ax = plot_edges(net['edge.conns'], net['node.coords'])
        # ax = plot_nodes(net['node.coords'], ax=ax)


if __name__ == '__main__':
    t = SKGROperationsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
