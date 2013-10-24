
"""
module Viscosity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(fluid,value=0.001,**params):
    return value

def na(fluid,**params):
    return 'n/a'

def Reynolds(fluid,uo=1,b=1,**params):
    T = fluid.pore_conditions['temperature']
    mu = uo*sp.exp(-1*b*T)
    return mu