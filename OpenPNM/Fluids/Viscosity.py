
"""
module Viscosity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(network,value=0.001,**params):
    return value

def na(network,**params):
    return 'n/a'

def Reynolds(network,uo=1,b=1,**params):
    T = network.pore_conditions['temperature']
    mu = uo*sp.exp(-1*b*T)
    return mu