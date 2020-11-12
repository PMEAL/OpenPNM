r"""
"""
import numpy as np
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def blank(target):
    pass


def invert(target, prop):
    r"""
    Inverts the given array

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    prop : string
        Dictionary key pointing the values to invert

    """
    return 1.0/target[prop]


def fraction(target, numerator, denominator):
    r"""
    Calculates the ratio between two values

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    numerator : string
        Dictionary key pointing the numerator values

    denominator : string
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
    target : OpenPNM Object
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
    target : OpenPNM Object
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
    target : OpenPNM Object
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
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    value : scalar
        The numerical value to apply

    Returns
    -------
    value : NumPy ndarray
        Array containing constant values equal to ``value``.

    Notes
    -----
    This model is mostly useless and for testing purposes, but might be used
    to 'reset' an array back to a default value.

    """
    return value


def product(target, prop1, prop2, **kwargs):
    r"""
    Calculates the product of multiple property values

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    prop1 : string
        The name of the first argument

    prop2 : string
        The name of the second argument

    Returns
    -------
    value : NumPy ndarray
        Array containing product values of ``target[prop1]``, ``target[prop2]``

    Notes
    -----
    Additional properties can be specified beyond just ``prop1`` and ``prop2``
    by including additional arguments in the function call (i.e. ``prop3 =
    'pore.foo'``).

    """
    value = target[prop1]*target[prop2]
    for item in kwargs.values():
        value *= target[item]
    return value


def scaled(target, prop, factor):
    r"""
    Scales an existing value by a factor.  Useful for constricting some throat
    property.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    prop : string
        The dictionary key of the array containing the values to be scaled.

    factor : scalar
        The factor by which the values should be scaled.

    Returns
    -------
    value : NumPy ndarray
        Array containing ``target[prop]`` values scaled by ``factor``.

    """
    value = target[prop]*factor
    return value
