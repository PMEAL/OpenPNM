"""
Basic Math
==========

"""
import logging
import numpy as np
logger = logging.getLogger(__name__)


__all__ = [
    'blank',
    'clip',
    'constant',
    'difference',
    'fraction',
    'invert',
    'normalize',
    'scaled',
    'summation',
    'product',
]


def blank(target):
    """Blank model used in PNM format, acting as a placeholder"""
    pass


def invert(target, prop):
    r"""
    Inverts the given array

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    prop : str
        Dictionary key pointing the values to invert

    """
    return 1.0/target[prop]


def difference(target, props):
    r"""
    Subtracts elements 1:N in `props` from element 0

    Parameters
    ----------
    target : OpenPNM dict
        The object to which the model is associated
    props : list
        A list of dict keys containing the values to operate on.  If the first
        element is A, and the next are B and C, then the results is A - B - C.
    """
    A = target[props[0]]
    for B in props[1:]:
        A = A - target[B]
    return A


def fraction(target, numerator, denominator):
    r"""
    Calculates the ratio between two values

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    numerator : str
        Dictionary key pointing the numerator values
    denominator : str
        Dictionary key pointing the denominator values

    """
    x = target[numerator]
    y = target[denominator]
    return x/y


def summation(target, props=[]):
    r"""
    Sums the values in the given arrays

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    props : list of dictionary keys
        The dictionary keys pointing the arrays whose elements should be summed

    """
    vals = np.zeros_like(target[props[0]])
    for item in props:
        vals += target[item]
    return vals


def normalize(target, prop, xmin=0, xmax=1):
    r"""
    Normalizes the given array between the supplied limits

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    xmin : float
        Lower limit of the re-scaled data
    xmax : float
        Upper limit of the re-scaled data

    """
    vals = target[prop]
    # Scale to 0 to 1
    vals = (vals - vals.min())/(vals.max() - vals.min())
    vals = vals*(xmax - xmin) + xmin
    return vals


def clip(target, prop, xmax, xmin=0):
    r"""
    Clips the given array within the supplied limits

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    xmin : float
        Values below this limit will be replaced with ``xmin``.
    xmax : float
        Values above this limit will be replaced with ``xmax``.

    """
    vals = np.clip(target[prop], xmin, xmax)
    return vals


def constant(target, value):
    r"""
    Places the given constant value into the target object

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    value : scalar
        The numerical value to apply

    Returns
    -------
    value : ndarray
        Array containing constant values equal to ``value``.

    Notes
    -----
    This model is mostly useless and for testing purposes, but might be used
    to 'reset' an array back to a default value.

    """
    return value


def product(target, props):
    r"""
    Calculates the product of multiple property values

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    props : list[str]
        The name of the arguments to be used for the product.

    Returns
    -------
    value : ndarray
        Array containing product values of ``target[props[0]]``,
        ``target[props[1]]``, etc.

    """
    value = np.ones_like(target[props[0]])
    for item in props:
        value *= target[item]
    return value


def scaled(target, prop, factor):
    r"""
    Scales an existing value by a factor.

    Useful for constricting some throat property.

    Parameters
    ----------
    target : Base
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    prop : str
        The dictionary key of the array containing the values to be scaled.
    factor : str
        The factor by which the values should be scaled.

    Returns
    -------
    value : ndarray
        Array containing ``target[prop]`` values scaled by ``factor``.

    """
    value = target[prop]*factor
    return value
