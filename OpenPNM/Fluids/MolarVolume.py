
"""
module Diffusivity
===============================================================================

"""

import scipy as _sp

def set_molar_volume(c):
    return {'molar_volume': c}

def ideal_gas_law(T,P):
    c = P/(8.314*T)
    return {'molar_volume': c}

