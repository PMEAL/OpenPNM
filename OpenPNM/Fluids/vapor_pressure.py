
"""
module vapor_pressure
===============================================================================

"""
import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def na(fluid,network,propname,**params):
    r"""
    Assigns nonsensical, but numerical value of -1.  
    This ensurse stability of other methods 
    but introduces the possibility of being misused.
    """
    value = -1
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def Antoine(fluid,network,propname,A=8.07131,B=1730.63,C=233.426,**params):
    r"""
    Uses Antoine equation to estimate vapor pressure of a pure component

    Parameters
    ----------
    A, B, C :  float, array_like
            Antoine vapor pressure constants for pure compounds

    """
    T = network.get_pore_data(phase=fluid,prop='temperature')
    value = (10**(A-B/(C+T-273.15)))*1.333e2
    network.set_pore_data(phase=fluid,prop=propname,data=value)