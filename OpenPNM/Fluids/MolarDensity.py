
"""
module MolarDensity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(network,value=40.89,**params):
    return value

def na(network,**params):
    return 'n/a'

def ideal_gas(fluid,**params):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas 

    """     
    R = 8.314    
    T = fluid.pore_conditions['temperature']
    P = fluid.pore_conditions['pressure']
    c = P/(R*T)
    return c

