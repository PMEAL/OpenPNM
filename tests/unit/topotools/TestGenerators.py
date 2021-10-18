import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from openpnm import topotools


class GeneratorTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def test_cubic(self):
        net = op.topotools.generators.cubic([3, 3, 3], 1e-3)
        assert net['vert.coords'].shape[0] == 27
        assert net['edge.conns'].shape[0] == 54

    def test_fcc(self):
        net = op.topotools.generators.fcc([3, 3, 3], 1e-3)
        assert net['vert.coords'].shape[0] == 172
        assert net['edge.conns'].shape[0] == 900

    def test_bcc(self):
        net = op.topotools.generators.bcc([3, 3, 3], 1e-3)
        assert net['vert.coords'].shape[0] == 91
        assert net['edge.conns'].shape[0] == 414

    def test_delaunay(self):
        np.random.seed(0)
        net, tri = op.topotools.generators.delaunay(points=20, shape=[1, 1, 1])
        assert net['vert.coords'].shape[0] == 20
        assert net['edge.conns'].shape[0] == 90

    def test_gabriel(self):
        np.random.seed(0)
        net = op.topotools.generators.gabriel(points=20, shape=[1, 1, 1])
        assert net['vert.coords'].shape[0] == 20
        assert net['edge.conns'].shape[0] == 51

    def test_voronoi(self):
        np.random.seed(0)
        net, vor = op.topotools.generators.voronoi(points=20, shape=[1, 1, 1])
        assert net['vert.coords'].shape[0] == 60
        assert net['edge.conns'].shape[0] == 109

    def test_voronoi_delaunay_dual(self):
        np.random.seed(0)
        f = op.topotools.generators.voronoi_delaunay_dual
        net, tri, vor = f(points=20, shape=[1, 1, 1])
        assert net['vert.coords'].shape[0] == 80
        assert net['edge.conns'].shape[0] == 457

    def test_cubic_template(self):
        im = np.ones([50, 50], dtype=bool)
        im[25:, ...] = False
        net = op.topotools.generators.cubic_template(template=im)
        assert net['vert.coords'].shape[0] == 1250
        assert net['edge.conns'].shape[0] == 2425


if __name__ == '__main__':

    t = GeneratorTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
