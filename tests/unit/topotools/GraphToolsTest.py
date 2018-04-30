import openpnm as op
import numpy as np
from numpy.testing import assert_approx_equal
from openpnm import topotools


class GraphToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_istriu(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istriu(am)
        am = net.create_adjacency_matrix(triu=True)
        assert topotools.istriu(am)
        am = am.T
        assert not topotools.istriu(am)

    def test_istril(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istril(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.istril(am)
        am = am.T
        assert topotools.istril(am)

    def test_istriangular(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert not topotools.istriangular(am)
        am = net.create_adjacency_matrix(triu=True)
        assert topotools.istriangular(am)
        am = am.T
        assert topotools.istriangular(am)

    def test_issymmetric(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix(triu=False)
        assert topotools.issymmetric(am)
        am = net.create_adjacency_matrix(triu=True)
        assert not topotools.issymmetric(am)
        am = am.T
        assert not topotools.issymmetric(am)

    def test_am_to_im(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        am = net.create_adjacency_matrix()
        im = net.create_incidence_matrix()
        im2 = topotools.am_to_im(am)
        assert im == im2


if __name__ == '__main__':

    t = GraphToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
