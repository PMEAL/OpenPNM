
"""
module VaporPressure
===============================================================================

"""
import OpenPNM
import scipy as sp

def set_as(fluid=None,Pv=3000):
    Pv = sp.array(Pv)
    fluid.update({'vapor_pressure': Pv})

def Antoine(fluid=None,T=298,A=8.07131,B=1730.63,C=233.426):
    r"""
    Takes Antoine coefficients and returns the vapor pressure
    """
    Pv = 10**(A-B/(C+T))
    OpenPNM.Fluids.VaporPresure(fluid,Pv)