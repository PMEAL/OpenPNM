r"""
===============================================================================
Submodule -- miscillaneous
===============================================================================

Models for applying basic phase properties

"""
import scipy as _sp


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
        property to be used in the correlation.

    """
    x = target[prop]
    value = 0.0
    for i in range(0, len(a)):
        value += a[i]*x**i
    return value
