import numpy as np
import openpnm as op


class VoronoiTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_voronoi_cube_num_points(self):
        net = op.network.Voronoi(points=30, shape=[1, 1, 1])
        assert net.Np > 30
        assert not np.all(net.coords[:, -1] == 0.0)

    def test_voronoi_square_num_points(self):
        net = op.network.Voronoi(points=30, shape=[1, 1, 0])
        assert net.Np > 30
        assert np.all(net.coords[:, -1] == 0.0)

    def test_voronoi_square_2D_points(self):
        np.random.rand(0)
        pts = np.random.rand(30, 2)
        net = op.network.Voronoi(points=pts, shape=[1, 1, 0])
        assert net.Np > 30
        assert np.all(net.coords[:, -1] == 0.0)

    def test_voronoi_square_3D_points(self):
        np.random.rand(0)
        pts = np.random.rand(30, 3)
        net = op.network.Voronoi(points=pts, shape=[1, 1, 0])
        assert net.Np > 30
        assert np.all(net.coords[:, -1] == 0.0)

    def test_voronoi_delaunay_dual(self):
        np.random.seed(0)
        dual = op.network.DelaunayVoronoiDual(points=100, shape=[1, 1, 1])
        assert dual.num_pores('delaunay') == 100
        assert dual.num_pores('voronoi') == 567

    def test_find_throat_facets(self):
        np.random.seed(0)
        dual = op.network.DelaunayVoronoiDual(points=10, shape=[1, 1, 1])
        f = dual.find_throat_facets(throats=[1, 5])
        assert np.all(f[0] == [48, 49, 50, 55, 57])
        assert np.all(f[1] == [48, 33, 30, 49])

    def test_find_pore_hulls(self):
        np.random.seed(0)
        dual = op.network.DelaunayVoronoiDual(points=10, shape=[1, 1, 1])
        f = dual.find_pore_hulls(pores=[0, 5])
        assert np.all(f[0] == [12, 14, 15, 19, 20, 21, 30,
                               33, 35, 48, 49, 50, 55, 57])
        assert np.all(f[1] == [36, 37, 38, 39, 40, 41, 42,
                               43, 51, 58, 60, 61])


if __name__ == '__main__':

    t = VoronoiTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
