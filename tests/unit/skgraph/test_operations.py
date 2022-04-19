import numpy as np
from openpnm._skgraph import generators as gen
from openpnm._skgraph.visualization import plot_edges, plot_nodes
import openpnm._skgraph.operations as ops


class SKGROperationsTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_join(self):
        g1 = gen.cubic([3, 3, 3])
        g2 = gen.cubic([3, 3, 3])
        g2['node.coords'] += np.array([0, 0, 3])
        g3 = ops.join(g1, g2)
        assert g3['edge.conns'].shape[0] == 108

        g1 = gen.cubic([3, 3, 3])
        g2 = gen.cubic([3, 3, 3])
        g2['node.coords'] += np.array([0, 0, 3])
        g3 = ops.join(g1, g2, L_max=1.1)
        assert g3['edge.conns'].shape[0] == 117

        g1 = gen.cubic([3, 3, 3])
        g2 = gen.cubic([3, 3, 3])
        g2['node.coords'] += np.array([0, 0, 3])
        g3 = ops.join(g1, g2, L_max=1.9)
        assert g3['edge.conns'].shape[0] == 157

    def test_add_nodes(self):
        g = gen.cubic([3, 3, 3])
        g = ops.add_nodes(g, [4, 4, 4])
        assert g['node.coords'].shape[0] == 28
        assert np.all(g['node.coords'][-1, :] == [4, 4, 4])
        g['node.float'] = np.ones(28, dtype=float)
        g['node.int'] = np.ones(28, dtype=int)
        g['node.bool'] = np.ones(28, dtype=bool)
        g = ops.add_nodes(g, [5, 5, 5])
        assert g['node.float'].shape[0] == 29
        assert np.isnan(g['node.float'][-1])
        assert g['node.int'][-1] == -2147483648
        assert g['node.bool'][-1] == False

    def test_add_edges(self):
        g = gen.cubic([3, 3, 3])
        g = ops.add_edges(g, [2, 4])
        assert g['edge.conns'].shape[0] == 55
        assert np.all(g['edge.conns'][-1, :] == [2, 4])
        g['edge.float'] = np.ones(55, dtype=float)
        g['edge.int'] = np.ones(55, dtype=int)
        g['edge.bool'] = np.ones(55, dtype=bool)
        g = ops.add_edges(g, [2, 4])
        assert g['edge.float'].shape[0] == 56
        assert np.isnan(g['edge.float'][-1])
        assert g['edge.int'][-1] == -2147483648
        assert g['edge.bool'][-1] == False

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
