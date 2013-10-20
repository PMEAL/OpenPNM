
"""
module Viscosity
===============================================================================

"""
import OpenPNM
import scipy as sp

def set_as(fluid=None,visc=0.001):
    visc = sp.array(visc)
    fluid.update({'viscosity':visc})
