r"""
===============================================================================
Submodule -- vapor_pressure
===============================================================================

Methods for predicing the vapor pressure of pure species

"""
import scipy as _sp

def Antoine(fluid,
            A,B,C,
            **kwargs):
    r"""
    Uses Antoine equation to estimate vapor pressure of a pure component

    Parameters
    ----------
    A, B, C :  float, array_like
            Antoine vapor pressure constants for pure compounds.  Note that
            these constants are traditionally reported such that they give
            vapor pressurein mmHg.  This function converts pressure to Pascals.

    """
    T = fluid['pore.temperature']
    value = (10**(A-B/(C+T-273.15)))*1.333e2
    return value
