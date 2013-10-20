
"""
module Viscosity
===============================================================================

"""
import OpenPNM
import scipy as sp

def set_as(fluid=None,visc=0.001):
    visc = sp.array(visc,ndmin=1)
    fluid.update({'viscosity':visc})
