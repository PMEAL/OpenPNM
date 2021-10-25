import pytest
import numpy as np
import openpnm as op


class CubicTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_connectivity(self):
        clist = [6, 14, 18, 20, 26]
        conns_dict = {
            6: np.array([[0, 1],
                         [2, 3],
                         [4, 5],
                         [6, 7],
                         [0, 2],
                         [1, 3],
                         [4, 6],
                         [5, 7],
                         [0, 4],
                         [1, 5],
                         [2, 6],
                         [3, 7]]),
            14: np.array([[0, 1],
                          [2, 3],
                          [4, 5],
                          [6, 7],
                          [0, 2],
                          [1, 3],
                          [4, 6],
                          [5, 7],
                          [0, 4],
                          [1, 5],
                          [2, 6],
                          [3, 7],
                          [0, 7],
                          [1, 6],
                          [2, 5],
                          [3, 4]]),
            18: np.array([[0, 1],
                          [2, 3],
                          [4, 5],
                          [6, 7],
                          [0, 2],
                          [1, 3],
                          [4, 6],
                          [5, 7],
                          [0, 4],
                          [1, 5],
                          [2, 6],
                          [3, 7],
                          [0, 3],
                          [4, 7],
                          [1, 2],
                          [5, 6],
                          [0, 5],
                          [2, 7],
                          [1, 4],
                          [3, 6],
                          [0, 6],
                          [1, 7],
                          [2, 4],
                          [3, 5]]),
            20: np.array([[0, 3],
                          [4, 7],
                          [1, 2],
                          [5, 6],
                          [0, 5],
                          [2, 7],
                          [1, 4],
                          [3, 6],
                          [0, 6],
                          [1, 7],
                          [2, 4],
                          [3, 5],
                          [0, 7],
                          [1, 6],
                          [2, 5],
                          [3, 4]]),
            26: np.array([[0, 1],
                          [2, 3],
                          [4, 5],
                          [6, 7],
                          [0, 2],
                          [1, 3],
                          [4, 6],
                          [5, 7],
                          [0, 4],
                          [1, 5],
                          [2, 6],
                          [3, 7],
                          [0, 7],
                          [1, 6],
                          [2, 5],
                          [3, 4],
                          [0, 3],
                          [4, 7],
                          [1, 2],
                          [5, 6],
                          [0, 5],
                          [2, 7],
                          [1, 4],
                          [3, 6],
                          [0, 6],
                          [1, 7],
                          [2, 4],
                          [3, 5]])
        }
        for x in clist:
            net = op.network.Cubic(shape=[2, 2, 2], connectivity=x)
            np.testing.assert_allclose(net.conns, conns_dict[x])

    def test_invalid_connectivity(self):
        invalid_clist = [8, 12, 45]
        for x in invalid_clist:
            with pytest.raises(Exception):
                _ = op.network.Cubic(shape=[3, 4, 5], connectivity=x)


if __name__ == '__main__':

    t = CubicTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
