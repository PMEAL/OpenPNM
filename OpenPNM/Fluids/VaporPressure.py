
"""
module VaporPressure
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(fluid,value=3000,**params):
    return value

def na(fluid,**params):
    return 'n/a'

def Antoine(fluid,A=8.07131,B=1730.63,C=233.426,**params):
    T = fluid.pore_conditions['temperature']
    Pv = 10**(A-B/(C+T))
    return Pv