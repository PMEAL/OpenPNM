import numpy as np
from numba import njit
from scipy.stats import rankdata


__all__ = [
    'qupc_initialize',
    'qupc_update',
    'qupc_compress',
    'qupc_reduce',
]


@njit
def qupc_initialize(size):
    return np.arange(size, dtype=np.int_)


@njit
def qupc_update(arr, ind, val):
    if ind == val:
        arr[ind] = val
    else:
        # Update array and do path compression simultaneously
        while arr[ind] != val:
            arr[ind] = arr[val]
            ind = val
            val = arr[val]
    return arr


def qupc_compress(arr):
    temp = rankdata(arr, method='dense')
    arr[:] = temp
    arr -= 1
    return arr


@njit
def qupc_reduce(arr):
    for i in range(len(arr)-1, 0, -1):
        arr[i] = arr[arr[i]]
    return arr


if __name__ == '__main__':
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
    print(a)
