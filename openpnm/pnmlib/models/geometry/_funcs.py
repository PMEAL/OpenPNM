import numpy as np
from pnmlib.core import get_data, count


__all__ = [
    "random_seeds",
    "product",
    "constant",
]


def random_seeds(
    target,
    seed=None,
    num_range=[0, 1],
):
    Np = count(target, 'network/pore')
    lo, hi = num_range
    vals = (np.random.rand(Np) - lo)/(hi - lo) + lo
    return vals


def product(
    target,
    props,
):
    vals = get_data(target, props[0])
    for i in range(1, len(props)):
        vals = np.multiply(vals, get_data(target, props[i]))
    return vals


def constant(
    target,
    value,
):
    return value
