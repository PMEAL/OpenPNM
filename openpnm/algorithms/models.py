import scipy as sp


def standard_kinetics(target, quantity, prefactor, exponent):
    r"""

    """
    X = target[quantity]
    A = target[prefactor]
    b = target[exponent]

    S1 = A*b*X**(b - 1)
    S2 = A*X**b*(1 - b)
    S = sp.vstack((S1, S2)).T
    return S
