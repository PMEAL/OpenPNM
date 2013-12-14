
"""
module molar_density
===============================================================================

"""
import scipy as sp

def constant(fluid,value=40.89,**params):
    fluid.pore_conditions['molar_density'] = value

def na(fluid,**params):
    fluid.pore_conditions['molar_density'] = -1

def ideal_gas(fluid,**params):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    R = 8.314
    T = fluid.pore_conditions['temperature']
    P = fluid.pore_conditions['pressure']
    c = P/(R*T)
    fluid.pore_conditions['molar_density'] = c

