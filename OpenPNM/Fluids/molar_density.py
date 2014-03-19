r"""
===============================================================================
Submodule -- molar_density
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

def ideal_gas(fluid,network,propname,**params):
    r"""
    Uses ideal gas equation to estimate molar density of a pure gas

    """
    R = sp.constants.R
    T = network.get_pore_data(phase=fluid,prop='temperature')
    P = network.get_pore_data(phase=fluid,prop='pressure')
    value = P/(R*T)
    network.set_pore_data(phase=fluid,prop=propname,data=value)