import numpy as np
from openpnm._skgraph.queries import (
    qupc_initialize,
    qupc_reduce,
    qupc_update,
    qupc_compress,
)


class SKGRQUPCTest:
    def setup_class(self):
        pass

    def test_basic(self):
        a = qupc_initialize(10)
        qupc_update(a, 4, 2)
        qupc_update(a, 7, 4)
        qupc_update(a, 9, 6)
        qupc_update(a, 6, 2)
        qupc_update(a, 5, 9)
        assert np.all(a == [0, 1, 2, 3, 2, 6, 2, 2, 8, 6])
        qupc_reduce(a)
        assert np.all(a == [0, 1, 2, 3, 2, 2, 2, 2, 8, 2])
        qupc_update(a, 9, 9)
        qupc_update(a, 0, 1)
        qupc_update(a, 8, 0)
        assert np.all(a == [1, 1, 2, 3, 2, 2, 2, 2, 1, 9])
        qupc_reduce(a)
        qupc_compress(a)
        assert np.all(a == [0, 0, 1, 2, 1, 1, 1, 1, 0, 3])


if __name__ == '__main__':
    t = SKGRQUPCTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
