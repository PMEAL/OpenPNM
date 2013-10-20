
"""
module MolarDensity
===============================================================================

"""
import OpenPNM
import scipy as sp

def set_as(fluid=None,c=40.89):
    c = sp.array(c)
    fluid.update({'molar_density': c})

def ideal_gas_law(fluid,T=298,P=101325):
    c = P/(8.314*T)
    OpenPNM.Fluids.MolarDensity.set_as(fluid,c)

