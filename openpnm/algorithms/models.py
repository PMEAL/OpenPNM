import scipy as sp


def standard_kinetics(target, quantity, prefactor, exponent, rate=True):
    r"""

    """
    X = target[quantity]
    A = target[prefactor]
    b = target[exponent]

    if rate:
        values = A*(X**b)
    else:
        S1 = A*b*(X**(b - 1))
        S2 = A*(1 - b)*(X**b)
        values = sp.vstack((S1, S2)).T
    return values
