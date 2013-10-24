
"""
module MolarDensity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(fluid,value=40.89,**params):
    return value

def na(fluid,**params):
    return 'n/a'

def ideal_gas(fluid,R=8.314,**params):
    T = fluid.pore_conditions['temperature']
    P = fluid.pore_conditions['pressure']
    c = P/(R*T)
    return c

