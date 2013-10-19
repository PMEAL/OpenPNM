
"""
module MolarDensity
===============================================================================

"""

import scipy as _sp

def set_as(c=40.89):
    return {'molar_density': c}

def ideal_gas_law(T,P):
    c = P/(8.314*T)
    return {'molar_density': c}

