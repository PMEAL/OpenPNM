import numpy as np
import openpnm as op
from numpy.testing import assert_array_almost_equal


class DelaunayGabrielTest:

    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_delaunay_square_with_reflect(self):
        np.random.seed(0)
        shape = [1, 1, 0]
        tri = op.network.Delaunay(points=30, shape=shape, reflect=False)
        assert op.topotools.isoutside(network=tri, shape=shape).sum() == 0
        tri = op.network.Delaunay(points=30, shape=shape, reflect=True)
        assert op.topotools.isoutside(network=tri, shape=shape).sum() == 0

    def test_delaunay_cube_with_trim_reflect(self):
        np.random.seed(0)
        shape = [1, 1, 1]
        tri = op.network.Delaunay(points=30, shape=shape, reflect=False)
        assert op.topotools.isoutside(network=tri, shape=shape).sum() == 0
        tri = op.network.Delaunay(points=30, shape=shape, reflect=True)
        assert op.topotools.isoutside(network=tri, shape=shape).sum() == 0

    def test_delaunay_square_with_2D_points(self):
        np.random.seed(0)
        pts = np.random.rand(50, 2)
        tri = op.network.Delaunay(points=pts, shape=[1, 1, 0])
        assert tri.coords.shape == (50, 3)
        assert np.all(tri.coords[:, :2] == pts)

    def test_delaunay_square_with_3D_points(self):
        np.random.seed(0)
        pts = np.random.rand(50, 3)
        tri = op.network.Delaunay(points=pts, shape=[1, 1, 0])
        assert tri.coords.shape == (50, 3)
        assert np.all(tri.coords[:, :2] == pts[:, :2])
        assert np.all(tri.coords[:, -1] != pts[:, -1])
        assert np.all(tri.coords[:, -1] == 0.0)

    def test_delaunay_square_with_num_points(self):
        np.random.seed(0)
        tri = op.network.Delaunay(points=30, shape=[1, 1, 0])
        assert op.topotools.dimensionality(network=tri).sum() == 2

    def test_delaunay_cube_with_points(self):
        np.random.seed(0)
        pts = np.random.rand(50, 3)
        tri = op.network.Delaunay(points=pts, shape=[1, 1, 1])
        assert tri.coords.shape == (50, 3)
        assert np.all(tri.coords == pts)

    def test_delaunay_cube_with_num_points(self):
        np.random.seed(0)
        tri = op.network.Delaunay(points=30, shape=[1, 1, 1])
        assert op.topotools.dimensionality(network=tri).sum() == 3

    def test_delaunay_disk_with_2D_points(self):
        np.random.seed(0)
        rqz = np.random.rand(50, 3)*np.array([1, 2*np.pi, 1])
        pts = np.vstack(op._skgraph.tools.cyl2cart(*rqz.T)).T
        tri = op.network.Delaunay(points=pts[:, :2], shape=[1, 0])
        assert tri.coords.shape == (50, 3)
        assert_array_almost_equal(tri.coords[:, :2], pts[:, :2], decimal=15)
        assert np.all(tri.coords[:, -1] != pts[:, -1])
        assert np.all(tri.coords[:, -1] == 0.0)

    def test_delaunay_disk_with_3D_points(self):
        np.random.seed(0)
        rqz = np.random.rand(50, 3)*np.array([1, 2*np.pi, 1])
        pts = np.vstack(op._skgraph.tools.cyl2cart(*rqz.T)).T
        tri = op.network.Delaunay(points=pts, shape=[1, 1])
        assert tri.coords.shape == (50, 3)
        assert_array_almost_equal(tri.coords, pts, decimal=15)

    def test_delaunay_disk_with_num_points(self):
        np.random.seed(0)
        tri = op.network.Delaunay(points=30, shape=[1, 0])
        assert op.topotools.dimensionality(network=tri).sum() == 2

    def test_delaunay_cylinder_with_points(self):
        np.random.seed(0)
        rqz = np.random.rand(50, 3)*np.array([1, 2*np.pi, 1])
        pts = np.vstack(op._skgraph.tools.cyl2cart(*rqz.T)).T
        tri = op.network.Delaunay(points=pts, shape=[1, 1])
        assert tri.coords.shape == (50, 3)
        assert_array_almost_equal(tri.coords, pts, decimal=15)

    def test_delaunay_sphere_with_num_points(self):
        np.random.seed(0)
        tri = op.network.Delaunay(points=30, shape=[1])
        assert op.topotools.dimensionality(network=tri).sum() == 3


if __name__ == '__main__':

    t = DelaunayGabrielTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
