
"""
module ThermalConductivity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(network,value=0.001,**params):
    return value

def na(network,**params):
    return 'n/a'

def ThermalConductivity(network,uo=1,b=1,**params):
    T = network.pore_conditions['temperature']
    mu = T #This is not implemented yet
    return mu

