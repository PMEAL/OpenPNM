
"""
module MolarDensity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(c=40.89,**params):
    return c

def na(**params):
    return 'n/a'

def ideal_gas(network,R=8.314,**params):
    T = network.pore_conditions['temperature']
    P = network.pore_conditions['pressure']
    c = P/(R*T)
    return c

