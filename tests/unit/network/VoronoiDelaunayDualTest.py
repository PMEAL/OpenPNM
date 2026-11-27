import numpy as np
import openpnm as op
import scipy as sp
import platform


ver = tuple([int(i) for i in sp.__version__.split('.')])


class VoronoiTest:

    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_dual_square_trim_and_reflect(self):
        np.random.seed(0)
        shape = [1, 1, 0]
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=False, reflect=False)
        assert op.topotools.isoutside(network=net, shape=shape).sum() > 0
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=True, reflect=False)
        assert op.topotools.isoutside(network=net, shape=shape).sum() == 0
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=False, reflect=True)
        assert op.topotools.isoutside(network=net, shape=shape).sum() > 0
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=True, reflect=True)
        assert op.topotools.isoutside(network=net, shape=shape).sum() == 0

    def test_dual_cube_trim_and_reflect(self):
        np.random.seed(0)
        shape = [1, 1, 1]
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=False, reflect=False)
        assert op.topotools.isoutside(network=net, shape=shape).sum() > 0
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=True, reflect=False)
        assert op.topotools.isoutside(network=net, shape=shape).sum() == 0
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=False, reflect=True)
        assert op.topotools.isoutside(network=net, shape=shape).sum() > 0
        net = op.network.DelaunayVoronoiDual(points=30, shape=shape,
                                             trim=True, reflect=True)
        assert op.topotools.isoutside(network=net, shape=shape).sum() == 0

    def test_dual_square_num_points(self):
        np.random.seed(0)
        net = op.network.DelaunayVoronoiDual(points=30, shape=[1, 1, 0])
        assert net.Np > 30
        assert np.all(net.coords[:, -1] == 0.0)
        assert np.all(op.topotools.dimensionality(net) == [True, True, False])

    def test_dual_cube_num_points(self):
        np.random.seed(0)
        net = op.network.DelaunayVoronoiDual(points=30, shape=[1, 1, 1])
        assert net.Np > 30
        assert not np.all(net.coords[:, -1] == 0.0)
        assert np.all(op.topotools.dimensionality(net) == [True, True, True])

    def test_dual_disk_num_points(self):
        np.random.seed(0)
        net = op.network.DelaunayVoronoiDual(points=30, shape=[1, 0])
        assert np.all(net.coords[:, -1] == 0.0)

    def test_dual_cylinder_num_points(self):
        np.random.seed(0)
        net = op.network.DelaunayVoronoiDual(points=30, shape=[1, 1])
        assert not np.all(net.coords[:, -1] == 0.0)
        assert np.all(op.topotools.dimensionality(net) == [True, True, True])

    def test_dual_sphere_num_points(self):
        np.random.seed(0)
        net = op.network.DelaunayVoronoiDual(points=30, shape=[1])
        assert not np.all(net.coords[:, -1] == 0.0)
        assert np.all(op.topotools.dimensionality(net) == [True, True, True])

    def test_dual_delaunay_dual_check_points(self):
        np.random.seed(0)
        dual = op.network.DelaunayVoronoiDual(points=100, shape=[1, 1, 1])
        assert dual.num_pores('delaunay') == 100
        assert dual.num_pores('voronoi') == 567

    def test_find_throat_facets(self):
        np.random.seed(0)
        dual = op.network.DelaunayVoronoiDual(points=10, shape=[1, 1, 1])
        f = dual.find_throat_facets(throats=[1, 5])
        if platform.system() != "Darwin":
            assert np.all(f[0] == [48, 49, 50, 55, 57])
            assert np.all(f[1] == [48, 33, 30, 49])
        else:
            assert np.all(f[0] == [35, 16, 50, 54, 31])
            assert np.all(f[1] == [32, 50, 51, 31])

    def test_find_pore_hulls(self):
        np.random.seed(0)
        dual = op.network.DelaunayVoronoiDual(points=10, shape=[1, 1, 1])
        f = dual.find_pore_hulls(pores=[0, 5])
        if platform.system() != "Darwin":
            assert np.all(f[0] == [12, 14, 15, 19, 20, 21, 30,
                                   33, 35, 48, 49, 50, 55, 57])
            assert np.all(f[1] == [36, 37, 38, 39, 40, 41, 42,
                                   43, 51, 58, 60, 61])
        else:
            assert np.all(f[0] == [12, 14, 15, 16, 20, 21, 22,
                                   31, 32, 35, 50, 51, 54, 56])
            assert np.all(f[1] == [38, 39, 40, 41, 42, 43, 44,
                                   45, 46, 47, 48, 49])

    def test_dual_relaxation_reflect(self):
        for shape, ndim in [([1, 1, 1], 3), ([1, 1, 0], 2)]:
            np.random.seed(0)
            net = op.network.DelaunayVoronoiDual(
                points=30,
                shape=shape,
                trim=True,
                reflect=True,
                f=0.5,
                relaxation=2,
            )
            assert net.num_pores('delaunay') == 30
            assert op.topotools.isoutside(network=net, shape=shape).sum() == 0
            surf_vor = net['pore.surface'] & net['pore.voronoi']
            assert surf_vor.sum() > 0
            coords = net['pore.coords'][surf_vor, :ndim]
            on_face = np.any(
                np.isclose(coords, 0.0) | np.isclose(coords, 1.0),
                axis=1,
            )
            assert np.all(on_face)


if __name__ == '__main__':

    t = VoronoiTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
