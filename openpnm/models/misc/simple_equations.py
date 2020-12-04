r""""
"""
import numpy as np
from openpnm.utils import logging
logger = logging.getLogger(__name__)


def generic_function(target, prop, func, **kwargs):
    r"""
    Runs an arbitrary function on the given data

    This allows users to place a customized calculation into the automatated
    model regeneration pipeline.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    prop : string
        The dictionary key containing the array to be operated on

    func : Numpy function
        A handle to the function to apply

    kwargs : keyward arguments
        All arguments required by the specific Numpy function

    Returns
    -------
    result : NumPy ndarray
        Array containing func(target[prop], **kwargs).

    Examples
    --------
    The following example shows how to use a Numpy function, but any function
    can be used, as long as it returns an array object:

    >>> import openpnm as op
    >>> import numpy as np
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> geo = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> geo['pore.rand'] = np.random.rand(geo.Np)
    >>> geo.add_model(propname='pore.cos',
    ...               model=op.models.misc.generic_function,
    ...               func=np.cos,
    ...               prop='pore.rand')
    """
    values = target[prop]
    result = func(values, **kwargs)
    if not isinstance(result, np.ndarray):
        logger.warning('Given function must return a Numpy array')
    return result


def linear(target, m, b, prop):
    r"""
    Calculates a property as a linear function of a given property

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    m, b : floats
        Slope and intercept of the linear corelation

    prop : string
        The dictionary key containing the independent variable or phase
        property to be used in the correlation.

    Returns
    -------
    value : NumPy ndarray
        Array containing ``m * target[prop] + b``.

    """
    x = target[prop]
    value = m*x + b
    return value


def polynomial(target, a, prop, **kwargs):
    r"""
    Calculates a property as a polynomial function of a given property

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    a : array_like
        A list containing the polynomial coefficients, where element 0 in the
        list corresponds to a0 and so on.  Note that no entries can be skipped
        so 0 coefficients must be sent as 0.

    prop : string
        The dictionary key containing the independent variable or phase
        property to be used in the polynomial.

    Returns
    -------
    value : NumPy ndarray
        Array containing ``Pn(target[prop])``, where ``Pn`` is nth order
        polynomial with coefficients stored in ``a``.

    """
    x = target[prop]
    value = 0.0
    for i in range(0, len(a)):
        value += a[i]*x**i
    return value
