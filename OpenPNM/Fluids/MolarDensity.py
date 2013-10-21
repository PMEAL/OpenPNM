
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

def ideal_gas(T=298,P=101325,R=8.314,**params):
    c = P/(R*T)
    return c

