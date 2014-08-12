r"""
===============================================================================
Submodule -- density
===============================================================================

"""
import scipy as _sp

def IABWS(fluid,coefficients,
          **kwargs):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    T = fluid['pore.temperature']
    P = fluid['pore.pressure']
    value = 'some_function'
    return value
