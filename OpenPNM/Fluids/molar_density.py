
"""
module molar_density
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

def constant(fluid,network,value,**params):
    r"""
    Assigns specified constant value
    """
    fluid.pore_conditions[propname] = value

def na(fluid,network,**params):
    value = -1
    fluid.pore_conditions[propname] = value

def ideal_gas(fluid,network,**params):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    R = 8.314
    T = fluid.pore_conditions['temperature']
    P = fluid.pore_conditions['pressure']
    value = P/(R*T)
    fluid.pore_conditions[propname] = value