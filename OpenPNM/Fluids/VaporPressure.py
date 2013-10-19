
"""
module VaporPressure
===============================================================================

"""

import scipy as _sp

def set_as(Pv=3000):
    return {'vapor_pressure': Pv}

def Antoine(A,B,C,T=298,P=101326):
    r"""
    Takes Antoine coefficients and returns the vapor pressure as a K-value (Pv/P)
    """
    Pv = 10**(A-B/(C+T))
    K = Pv/P
    return {'vapor_pressure': K}