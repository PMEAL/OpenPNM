import numpy as np
from openpnm._skgraph import generators as gen
from openpnm._skgraph import settings
from openpnm._skgraph.visualization import plot_edges, plot_nodes
settings.node_prefix = 'node'
settings.edge_prefix = 'edge'


class SKGRGeneratorsTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_cubic(self):
        net = gen.cubic(shape=[3, 2, 1])
        assert net['node.coords'].shape[0] == 6
        assert net['edge.conns'].shape[0] == 7
        net = gen.cubic([3, 3, 3], 1e-3)
        assert net['node.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 54

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

    def test_fcc(self):
        net = gen.fcc([3, 3, 3], 1e-3, mode='triangulation')
        assert net['node.coords'].shape[0] == 63
        assert net['edge.conns'].shape[0] == 294
        net = gen.fcc([3, 3, 3], 1e-3, mode='kdtree')
        assert net['node.coords'].shape[0] == 63
        assert net['edge.conns'].shape[0] == 294

    def test_bcc(self):
        net = gen.bcc([3, 3, 3], 1e-3, mode='triangulation')
        assert net['node.coords'].shape[0] == 35
        assert net['edge.conns'].shape[0] == 130
        net = gen.bcc([3, 3, 3], 1e-3, mode='kdtree')
        assert net['node.coords'].shape[0] == 35
        assert net['edge.conns'].shape[0] == 130

    def test_delaunay(self):
        np.random.seed(0)
        net, tri = gen.delaunay(points=20, shape=[1, 1, 1])
        assert net['node.coords'].shape[0] == 20
        assert net['edge.conns'].shape[0] == 95

    def test_gabriel(self):
        np.random.seed(0)
        net = gen.gabriel(points=20, shape=[1, 1, 1])
        assert net['node.coords'].shape[0] == 20
        assert net['edge.conns'].shape[0] == 46

    def test_voronoi_cubic(self):
        np.random.seed(0)
        net, vor = gen.voronoi(points=20, shape=[1, 1, 1], trim=False)
        assert net['node.coords'].shape[0] == 65
        assert net['edge.conns'].shape[0] == 119
        np.random.seed(0)
        net, vor = gen.voronoi(points=20, shape=[1, 1, 1], trim=True)
        assert net['node.coords'].shape[0] == 39
        assert net['edge.conns'].shape[0] == 61

    def test_voronoi_square(self):
        np.random.seed(0)
        points = gen.tools.generate_base_points(50, [1, 1, 0], reflect=True)
        assert len(points) == 250
        net, vor = gen.voronoi(points=points, shape=[1, 1, 0], trim=True)
        assert len(net['node.coords']) == 100
        assert len(net['edge.conns']) == 146
        ax = plot_edges(net['edge.conns'], net['node.coords'])
        ax = plot_nodes(net['node.coords'], ax=ax)

    def test_voronoi_delaunay_dual(self):
        np.random.seed(0)
        net, tri, vor = gen.voronoi_delaunay_dual(points=20,
                                                  shape=[1, 1, 1],
                                                  trim=False)
        assert net['node.coords'].shape[0] == 85
        assert net['edge.conns'].shape[0] == 491
        np.random.seed(0)
        net, tri, vor = gen.voronoi_delaunay_dual(points=20,
                                                  shape=[1, 1, 1],
                                                  trim=True)
        assert net['node.coords'].shape[0] == 59
        assert net['edge.conns'].shape[0] == 115

    def test_cubic_template(self):
        im = np.ones([50, 50], dtype=bool)
        im[25:, ...] = False
        net = gen.cubic_template(template=im)
        assert net['node.coords'].shape[0] == 1250
        assert net['edge.conns'].shape[0] == 2425


if __name__ == '__main__':
    t = SKGRGeneratorsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
