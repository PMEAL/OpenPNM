r"""
===============================================================================
Submodule -- molar_density
===============================================================================

"""
import scipy as sp

def ideal_gas(phase,**kwargs):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    R = sp.constants.R
    T = phase['pore.temperature']
    P = phase['pore.pressure']
    value = P/(R*T)
    return value
