import openpnm as op
import scipy as sp


class DelaunayGabrielTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_gabriel_and_delaunay_cubic(self):
        sp.random.seed(0)
        dn = op.network.Delaunay(shape=[1, 1, 1], num_points=50)
        sp.random.seed(0)
        gn = op.network.Gabriel(shape=[1, 1, 1], num_points=50)
        assert gn.Nt < dn.Nt
        assert gn.num_pores('internal') == 50
        assert dn.num_pores('internal') == 50
        assert gn.num_pores('surface') == 75

    def test_gabriel_and_delaunay_square(self):
        sp.random.seed(0)
        dn = op.network.Delaunay(shape=[1, 1, 0], num_points=50)
        sp.random.seed(0)
        gn = op.network.Gabriel(shape=[1, 1, 0], num_points=50)
        assert gn.Nt < dn.Nt
        assert gn.num_pores('internal') == 50
        assert dn.num_pores('internal') == 50
        assert gn.num_pores('surface') == 24


if __name__ == '__main__':

    t = DelaunayGabrielTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
