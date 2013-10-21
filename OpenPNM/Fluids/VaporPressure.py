
"""
module VaporPressure
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(network,Pv=3000,**params):
    return Pv

def na(network,**params):
    return 'n/a'

def Antoine(network,A=8.07131,B=1730.63,C=233.426,**params):
    T = network.pore_conditions['temperature']
    Pv = 10**(A-B/(C+T))
    return Pv