r"""
===============================================================================
Submodule -- miscillaneous
===============================================================================

Models for applying basic phase properties

"""
import scipy as sp

def constant(phase,value,**kwargs):
    r"""
    Assigns specified constant value
    """
    temp = sp.ones(sp.shape(phase.pores()))*value
    return temp

def random(phase,seed=None,**kwargs):
    r"""
    Assigns specified constant value
    """
    sp.random.seed(seed)
    value = sp.random.rand(sp.shape(phase.pores())[0])
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
    
    
    
    
    
    
    
    
    