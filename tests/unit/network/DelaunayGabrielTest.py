import numpy as np
import scipy as sp
import openpnm as op


class DelaunayGabrielTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_gabriel_and_delaunay_cubic(self):
        np.random.seed(0)
        dn = op.network.Delaunay(shape=[1, 1, 1], points=50, trim=True)
        np.random.seed(0)
        gn = op.network.Gabriel(shape=[1, 1, 1], points=50, trim=True)
        assert gn.Nt < dn.Nt
        assert gn.num_pores(['internal', 'surface'], mode='union') == 50
        assert dn.num_pores(['internal', 'surface'], mode='union') == 50
        assert gn.num_pores('boundary') == 75

    def test_gabriel_and_delaunay_square(self):
        np.random.seed(0)
        dn = op.network.Delaunay(shape=[1, 1, 0], points=50)
        np.random.seed(0)
        gn = op.network.Gabriel(shape=[1, 1, 0], points=50)
        assert gn.Nt < dn.Nt
        assert gn.num_pores(['internal', 'surface'], mode='union') == 50
        assert dn.num_pores(['internal', 'surface'], mode='union') == 50
        assert gn.num_pores('boundary') == 24

    def test_add_boundary_pores(self):
        np.random.seed(0)
        dn = op.network.Delaunay(shape=[1, 1, 1], points=50)
        dn.add_boundary_pores(offset=0.1)
        assert np.all(np.amin(dn['pore.coords']) == -0.1)
        assert np.all(np.amax(dn['pore.coords']) == 1.1)


if __name__ == '__main__':

    t = DelaunayGabrielTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
