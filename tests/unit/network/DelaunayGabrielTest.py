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
        dn = op.network.Delaunay(shape=[1, 1, 1], points=50)
        # assert dn.num_pores(['internal', 'surface'], mode='union') == 50

    def test_gabriel_and_delaunay_square(self):
        np.random.seed(0)
        dn = op.network.Delaunay(shape=[1, 1, 0], points=50)
        # assert dn.num_pores(['internal', 'surface'], mode='union') == 50


if __name__ == '__main__':

    t = DelaunayGabrielTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
