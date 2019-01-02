r"""

.. autofunction:: openpnm.models.misc.basic_math.constant
.. autofunction:: openpnm.models.misc.basic_math.product
.. autofunction:: openpnm.models.misc.basic_math.scaled

"""

import numpy as np
import scipy.stats as spts
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def constant(target, value):
    r"""
    Places a constant value into the target object

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    value : scalar
        The numerical value to apply

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
    """
    value = target[prop]*factor
    return value
