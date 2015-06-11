r"""
===============================================================================
Submodule -- miscillaneous
===============================================================================

Models for applying basic phase properties

"""
import scipy as _sp


def constant(phase, value, **kwargs):
    r"""
    Assigns specified constant value
    """
    temp = _sp.ones(_sp.shape(phase.pores()))*value
    return temp


def random(phase, seed=None, **kwargs):
    r"""
    Assigns specified constant value
    """
    _sp.random.seed(seed=seed)
    value = _sp.random.rand(phase.Np)
    return value


def linear(phase, m, b, poreprop, **kwargs):
    r"""
    Calculates a property as a linear function of a given property

    Parameters
    ----------
    m, b : floats
        Slope and intercept of the linear corelation

    poreprop : string
        The dictionary key containing the independent variable or phase
        property to be used in the correlation.

    """
    T = phase[poreprop]
    value = b + m*T
    return value


def polynomial(phase, a, poreprop, **kwargs):
    r"""
    Calculates a property as a polynomial function of a given property

    Parameters
    ----------
    a : array_like
        A list containing the polynomial coefficients, where element 0 in the
        list corresponds to a0 and so on.  Note that no entries can be skipped
        so 0 coefficients must be sent as 0.

    poreprop : string
        The dictionary key containing the independent variable or phase
        property to be used in the correlation.

    """
    x = phase[poreprop]
    value = 0.0
    for i in range(0, len(a)):
        value += a[i]*x**i
    return value


def ideal_mixture(phase, poreprop,
                  composition='pore.mole_fraction',
                  **kwargs):
    r"""
    Calcualtes a given mixture property as the composition weighted average
    of the pure compononent properties

    Parameters
    ----------
    poreprop : string
        The dictionary key containing the independent variable or phase
        property to be used in the correlation.

    composition : string, optional (default is 'pore.mole_fraction')
        The name of the pore property where the composition information
        is stored on each pure component

    Returns
    -------
    The composition weighted average of the given property

    Notes
    -----
    The average is calculated as follows:

    .. math::

        P_{mixture}=\Sigma(x_{i}P_{i})

    where

        :math:`P_{mix}` is the average mixture property

        :math:`x_{i}` is the fractional composition of species *i*

        :math:`P_{i}` is the property of interest for pure species *i*


    """
    value = _sp.zeros((phase.Np,))
    for comp in phase._phases:
        value = value + comp[poreprop]*comp[composition]
    return value


def mixture_value(phase, propname, **kwargs):
    r"""
    Adopts the specified property value from the parent mixture phase

    Parameters
    ----------
    propname :
        The propname to which this model is assigned (i.e. 'pore.temperature')
        is automatically passed and used as the property name to fetch from
        the mixture object

    """
    mixture = phase._phases[0]
    vals = mixture[propname]
    return vals
