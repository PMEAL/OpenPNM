r"""
===============================================================================
Submodule -- molar_density
===============================================================================

"""
import scipy as _sp

def ideal_gas(fluid,**kwargs):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    R = _sp.constants.R
    T = fluid['pore.temperature']
    P = fluid['pore.pressure']
    value = P/(R*T)
    return value
