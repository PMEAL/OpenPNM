r"""
===============================================================================
Submodule -- miscillaneous
===============================================================================

Models for applying basic phase properties

"""
import scipy as _sp

def constant(phase,value,**kwargs):
    r"""
    Assigns specified constant value
    """
    temp = _sp.ones(_sp.shape(phase.pores()))*value
    return temp

def random(phase,seed=None,**kwargs):
    r"""
    Assigns specified constant value
    """
    _sp.random.seed(seed)
    value = _sp.random.rand(_sp.shape(phase.pores())[0])
    return value

def linear(phase,m,b,poreprop='pore.temperature',**kwargs):
    r"""
    Calculates a property as a linear function of a given property

    Parameters
    ----------
    m, b: floats
        Slope and intercept of the linear corelation

    poreprop : string
        The property name of the independent variable or phase property.  The
        default is 'pore.temperature'.

    """
    T = phase[poreprop]
    value = b + m*T
    return value

def polynomial(phase,a,poreprop='pore.temperature',**kwargs):
    r"""
    Calculates a property as a polynomial function of a given property

    Parameters
    ----------
    a: list of floats
        A list containing the polynomial coefficients, where element 0 in the
        list corresponds to a0 and so on.  Note that no entries can be skipped
        so 0 coefficients must be sent as 0.

    poreprop : string
        The property name of the independent variable or phase property.  The
        default is 'pore.temperature'.

    """
    x = phase[poreprop]
    value = 0.0
    for i in range(0,len(a)):
        value += a[i]*x**i
    return value

def ideal_mixture(phase,propname,pore_prop=None,composition='pore.mole_fraction',**kwargs):
    r'''
    Calcualtes a given mixture property as the composition weighted average
    of the pure compononent properties

    Parameters
    ----------
    pore_prop : string
        The name of the property on the pure component
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


    '''
    if pore_prop == None:
        pore_prop = propname
    value = _sp.zeros((phase.Np,))
    for comp in phase._phases:
        value= value + comp[pore_prop]*comp[composition]
    return value





