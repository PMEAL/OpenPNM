import numpy as np
import openpnm as op


class VoronoiTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_voronoi_num_points(self):
        net = op.network.Voronoi(points=30, shape=[1, 1, 1])
        assert net.Np > 30

    def test_voronoi_points(self):
        points = np.random.rand(30, 3)
        net = op.network.Voronoi(points=points, shape=[1, 1, 1])
        assert net.Np > 30

    def test_voronoi_2d(self):
        net = op.network.Voronoi(points=30, shape=[1, 1, 0])
        assert net.Np > 30


if __name__ == '__main__':

    t = VoronoiTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
