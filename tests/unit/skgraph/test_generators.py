import numpy as np
from openpnm._skgraph import generators as gen
from openpnm._skgraph import operations as ops
from openpnm._skgraph.operations import trim_nodes
from openpnm._skgraph.tools import isoutside


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
        net, tri = gen.delaunay(points=20, shape=[1, 1, 1], reflect=False)
        assert net['node.coords'].shape[0] == 20
        assert net['edge.conns'].shape[0] == 95

    def test_gabriel(self):
        np.random.seed(0)
        net = gen.gabriel(points=20, shape=[1, 1, 1], reflect=False)
        assert net['node.coords'].shape[0] == 20
        assert net['edge.conns'].shape[0] == 46

    def test_voronoi_cubic(self):
        np.random.seed(0)
        net, vor = gen.voronoi(points=20, shape=[1, 1, 1], reflect=False, trim=False)
        assert net['node.coords'].shape[0] == 65
        assert net['edge.conns'].shape[0] == 119
        np.random.seed(0)
        net, vor = gen.voronoi(points=20, shape=[1, 1, 1], reflect=False, trim=True)
        assert net['node.coords'].shape[0] == 39
        assert net['edge.conns'].shape[0] == 61

    def test_voronoi_square(self):
        np.random.seed(0)
        points = gen.tools.generate_base_points(50, [1, 1, 0], reflect=True)
        assert len(points) == 250
        net, vor = gen.voronoi(points=points, shape=[1, 1, 0], trim=True)
        assert net['node.coords'].shape[0] == 102
        assert net['edge.conns'].shape[0] == 151

    def test_voronoi_circle(self):
        np.random.seed(0)
        shape = [4, 0]
        points = gen.tools.generate_base_points(1000, shape, reflect=True)
        assert len(points) == 2000
        net, vor = gen.voronoi(points=points, shape=shape, trim=False)
        trim = isoutside(net, shape=shape, rtol=0.0015)
        net = trim_nodes(net, np.where(trim)[0])
        assert net['node.coords'].shape[0] == 2098
        assert net['edge.conns'].shape[0] == 3097

    def test_voronoi_cylinder(self):
        np.random.seed(0)
        shape = [2, 4]
        points = gen.tools.generate_base_points(700, shape, reflect=True)
        net, vor = gen.voronoi(points=points, shape=shape, trim=False)
        trim = isoutside(net, shape=shape, rtol=[0.05, 0])
        net = trim_nodes(net, np.where(trim)[0])
        assert net['node.coords'].shape[0] == 5800
        assert net['edge.conns'].shape[0] == 11272

    def test_voronoi_sphere(self):
        np.random.seed(0)
        shape = [2]
        points = gen.tools.generate_base_points(700, shape, reflect=True)
        net, vor = gen.voronoi(points=points, shape=shape, trim=False)
        trim = isoutside(net, shape=shape, rtol=[0.05])
        net = trim_nodes(net, np.where(trim)[0])
        assert net['node.coords'].shape[0] == 5656
        assert net['edge.conns'].shape[0] == 10953

    def test_voronoi_delaunay_dual_square(self):
        np.random.seed(0)
        shape = [1, 1, 0]
        points = gen.tools.generate_base_points(50, shape, reflect=True)
        assert len(points) == 250
        net, vor, tri = gen.voronoi_delaunay_dual(points=points,
                                                  shape=shape,
                                                  reflect=False,
                                                  trim=True)
        assert net['node.coords'].shape[0] == 152
        assert net['edge.conns'].shape[0] == 552
        net = ops.trim_nodes(net, net['node.delaunay'])
        assert net['node.coords'].shape[0] == 102
        assert net['edge.conns'].shape[0] == 151

    def test_voronoi_delaunay_dual(self):
        np.random.seed(0)
        net, vor, tri = gen.voronoi_delaunay_dual(points=20,
                                                  shape=[1, 1, 1],
                                                  trim=False,
                                                  reflect=False)
        assert net['node.coords'].shape[0] == 85
        assert net['edge.conns'].shape[0] == 491
        np.random.seed(0)
        net, vor, tri = gen.voronoi_delaunay_dual(points=20,
                                                  shape=[1, 1, 1],
                                                  trim=True,
                                                  reflect=False)
        assert net['node.coords'].shape[0] == 59
        assert net['edge.conns'].shape[0] == 311
        np.random.seed(0)
        net, vor, tri = gen.voronoi_delaunay_dual(points=20,
                                                  shape=[1, 1, 1],
                                                  trim=True,
                                                  reflect=True)
        assert net['node.coords'].shape[0] == 124
        assert net['edge.conns'].shape[0] == 597

    def test_cubic_template(self):
        im = np.ones([50, 50], dtype=bool)
        im[25:30, 25:40] = False
        net = gen.cubic_template(template=im)
        assert net['node.coords'].shape[0] == 2425
        assert net['edge.conns'].shape[0] == 4730


# ax = plot_edges(net)
# ax = plot_nodes(net, ax=ax)


if __name__ == '__main__':
    t = SKGRGeneratorsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
