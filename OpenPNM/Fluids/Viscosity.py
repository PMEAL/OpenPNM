
"""
module Viscosity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(visc=0.001,**params):
    return visc

def na(**params):
    return 'n/a'

def Reynolds(network,uo=1,b=1,**params):
    T = network.pore_conditions['temperature']
    visc = uo*sp.exp(-1*b*T)
    return visc