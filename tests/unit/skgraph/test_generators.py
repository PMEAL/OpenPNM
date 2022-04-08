from openpnm._skgraph import generators as gen
from openpnm._skgraph import settings
settings.node_prefix = 'node'
settings.edge_prefix = 'edge'


class CubicTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_cubic(self):
        net = gen.cubic(shape=[3, 2, 1])
        assert net['node.coords'].shape[0] == 6
        assert net['edge.conns'].shape[0] == 7

    def test_cubic_w_connectivity(self):
        shape = [3, 3, 3]
        net = gen.cubic(shape=shape, connectivity=14)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 86
        net = gen.cubic(shape=shape, connectivity=18)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 126
        net = gen.cubic(shape=shape, connectivity=20)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 104
        net = gen.cubic(shape=shape, connectivity=26)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 158


if __name__ == '__main__':
    t = CubicTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
