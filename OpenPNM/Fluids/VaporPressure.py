
"""
module VaporPressure
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(Pv=3000,**params):
    return Pv

def na(**params):
    return 'n/a'

def Antoine(T=298,A=8.07131,B=1730.63,C=233.426,**params):
    Pv = 10**(A-B/(C+T))
    return Pv